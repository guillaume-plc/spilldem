/***************************************************************
#                         spillDEM v0.1                        #
****************************************************************
#                                                              #
#     Fast DEM surface depressions filling using the spill     #
#   elevation and least-cost search approach from              #
#   [Wang, Lei & Liu, Holiday. (2006)]                         #
#   (http://dx.doi.org/10.1080/13658810500433453).             #
#                                                              #
#     Basically a standalone version of SAGA's implementation, #
#   relying only on GDAL for IO (see [SAGA's source code]      #
#   (https://sourceforge.net/projects/saga-gis/) for the       #
#   'Fill sinks (Wang & Liu)' preprocessing algorithm). As     #
#   for SAGA's implementation the original algorithm is        #
#   modified to allow for preservation of a minimum slope      #
#   gradient between cells.                                    #
#                                                              #
***************************************************************/

#include <getopt.h>
#include <iostream>
#include <queue>
#include <climits>
#include <cmath>
#include "gdal_priv.h"
#include "cpl_conv.h"

#include "SpillDEM.h" // config file

static void usage(const char* name)
{
    printf("%s version %d.%d\n"
           "usage: %s <options> datasource\n"
           "Options:\n"
            "\t-o, --output        filled DEM output file\n"
            "\t-f, --flow          D8 flow direction output file\n"
            "\t-m, --minslope      minimum preserved slope gradient"
            "\t-v, --verbose       display information messages\n"
            "\n"
            "\t-h, --help          display this message and exit\n",
            name, SpillDEM_VERSION_MAJOR, SpillDEM_VERSION_MINOR, name);
}

struct node
{
    float spill;
    int x;
    int y;

    node(float spill, int x, int y) 
        : spill(spill), x(x), y(y)
    {}

    bool operator<(const node& rhs) const
    {
        return spill > rhs.spill;
    }
};

struct dir
{
    int dx;
    int dy;

    dir(int dx, int dy)
        : dx(dx), dy(dy)
    {}
};

int main(int argc, char* argv[])
{
    const option long_opts[] =
    {
        {"output", required_argument, nullptr, 'o'},
        {"flow", required_argument, nullptr, 'f'},
        {"minslope", required_argument, nullptr, 'm'},
        {"verbose", no_argument, nullptr, 'v'},
        {"help", no_argument, nullptr, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    bool verbose = false;
    float minslope = 0.1;
    std::string infile = "";
    std::string spill_outfile = "filled.tif";
    std::string flow_outfile = "flow.tif";
    while ((opt = getopt_long(argc, argv, ":o:f:m:vh", long_opts, nullptr)) != -1) 
    {
        switch (opt) 
        {
        case 'v':
            verbose = true;
            break;
        case 'o':
            spill_outfile = std::string(optarg);
            break;
        case 'f':
            flow_outfile = std::string(optarg);
            break;
        case 'm':
            minslope = std::atof(optarg);
            break;
        case 'h':
            usage(argv[0]);
            exit(EXIT_SUCCESS);
            break;
        case '?':
            usage(argv[0]);
            fprintf(stderr, "Error: Unknown option -%c\n", (char)optopt);
            exit(EXIT_FAILURE);
            break;
        case ':':
            usage(argv[0]);
            fprintf(stderr, "Error: Option -%c requires an argument\n", (char)optopt);
            exit(EXIT_FAILURE);
            break;
        }
    }

    if (optind >= argc)
    {
        usage(argv[0]);
        fprintf(stderr, "Error: No data source specified.\n");
        exit(EXIT_FAILURE);
    }
    infile = argv[optind];

    GDALAllRegister();
    const char *format = "GTiff";
    GDALDriver *driver;
    driver = GetGDALDriverManager()->GetDriverByName(format);
    if( driver == nullptr )
    {
        exit(EXIT_FAILURE);
    }

    GDALDataset *srcDataset;
    srcDataset = (GDALDataset *)GDALOpen(infile.c_str(), GA_ReadOnly);
    if ( srcDataset == nullptr)
    {
        exit(EXIT_FAILURE);
    }

    GDALRasterBand *srcBand, *flowBand, *spillBand;
    srcBand = srcDataset->GetRasterBand(1);
    const int xSize = srcBand->GetXSize(), ySize = srcBand->GetYSize();
    double nodata = srcBand->GetNoDataValue();

    GDALDataset *flowDataset, *spillDataset;
    flowDataset = driver->Create(flow_outfile.c_str(), xSize, ySize, 1, GDT_Byte, NULL);
    spillDataset = driver->Create(spill_outfile.c_str(), xSize, ySize, 1, GDT_Float32, NULL);
    if ( flowDataset == nullptr )
    {
        GDALClose(srcDataset);
        exit(EXIT_FAILURE);
    }
    if ( spillDataset == nullptr )
    {
        GDALClose(srcDataset);
        GDALClose(flowDataset);
        exit(EXIT_FAILURE);
    }
    double adfGeoTransform[6];
    srcDataset->GetGeoTransform(adfGeoTransform);
    flowDataset->SetGeoTransform(adfGeoTransform);
    flowDataset->SetSpatialRef(srcDataset->GetSpatialRef());
    spillDataset->SetGeoTransform(adfGeoTransform);
    spillDataset->SetSpatialRef(srcDataset->GetSpatialRef());
    
    flowBand = flowDataset->GetRasterBand(1);
    flowBand->SetNoDataValue(255);
    spillBand = spillDataset->GetRasterBand(1);
    spillBand->SetNoDataValue(nodata);

    float *elev;
    elev = (float *) CPLMalloc(sizeof(float)*xSize*ySize);
    srcBand->RasterIO(GF_Read, 0, 0, xSize, ySize, elev, xSize, ySize, GDT_Float32, 0, 0);

    bool preserve;
    std::array<dir, 8> ngh = { dir(1, 0), dir(1, -1), dir(0, -1),
                               dir(-1, -1), dir(-1, 0), dir(-1, 1),
                               dir(0, 1), dir(1, 1) };
    float pixelSizeX = adfGeoTransform[1], pixelSizeY = adfGeoTransform[5];
    float diaglength = std::sqrt(pixelSizeX * pixelSizeX + pixelSizeY * pixelSizeY);
    std::array<float, 8> length = { pixelSizeX, diaglength, pixelSizeY,
                                    diaglength, pixelSizeX, diaglength,
                                    pixelSizeY, diaglength};
    std::array<float, 8> mindiff;

    if( minslope > 0.0 )
	{
		minslope = std::tan(minslope * M_PI / 180.0);
		for(int d=0; d<8; d++)
			mindiff[d] = minslope * length[d];
		preserve = true;
	}
	else
    {
		preserve = false;
    }

    auto getNeighbourX = [&](int x, int d){ return x + ngh[d].dx; };
    auto getNeighbourY = [&](int y, int d){ return y + ngh[d].dy; };
    auto getIndex = [&](int x, int y){ return y * xSize + x; };
    auto isInBounds = [&](int x, int y){ return x >= 0 && x < xSize && y >= 0 && y < ySize; };

    std::array<unsigned char, 9> d8 = {1, 2, 4, 8, 16, 32, 64, 128, 0};
    std::array<unsigned char, 9> ldd = {6, 3, 2, 1, 4, 7, 8, 9, 0};
    std::priority_queue<node> queue;
    std::vector<bool> queued(xSize*ySize, false);
    std::vector<bool> processed(xSize*ySize, false);
    std::vector<char> flowdir(xSize*ySize, 0);

    auto getFlowDir = [&](int x, int y, int z)
    {
        float maxgrad = -1.0, grad;
        char dmax = 8;
        int nx, ny, n;
        for (int d = 0; d < 8; d++)
        {
            nx = getNeighbourX(x, d);
            ny = getNeighbourY(y, d);
            n = getIndex(nx ,ny);
            if ( isInBounds(nx, ny) && processed[n] && elev[n] <= z)
            {
                grad = (z - elev[n]) / length[d];
                if (grad > maxgrad)
                {
                    maxgrad = grad;
                    dmax = d;
                }
            }
        }
        return dmax;
    };

    int c, n, nx, ny;
    float z, nz;

    // Initialize edge cells
    for (int x = 0; x < xSize; x++)
    {
        for (int y = 0; y < ySize; y++)
        {
            int n = getIndex(x, y);
            z = elev[n];
            if (elev[n] == nodata)
            {
                processed[n] = true;
                flowdir[n] = 255;
            }
            else
            {
                for (int d = 0; d < 8; d++)
                {
                    nx = getNeighbourX(x, d);
                    ny = getNeighbourY(y, d);
                    if ( !isInBounds(nx, ny) || elev[getIndex(nx, ny)] == nodata )
                    {
                        flowdir[n] = 255;
                        queue.push(std::move(node(z, x, y)));
                        queued[n] = true;
                        break;
                    }
                }
            }
        }
    }

    node current(0.0f, 0, 0);
    float maxgrad, grad;
    char dmax;
    while (!queue.empty())
    {
        current = queue.top();
        queue.pop();
        c = getIndex(current.x, current.y);
        processed[c] = true;
        queued[c] = false;
        maxgrad = -1.0;
        dmax = 8;
        z = current.spill;
        for (int d = 0; d < 8; d++)
        {
            nx = getNeighbourX(current.x, d);
            ny = getNeighbourY(current.y, d);
            n = getIndex(nx, ny);
            if ( isInBounds(nx, ny) && !queued[n])
            {
                /*
                if (processed[n])
                {
                    if (!flowdir[c] && (elev[c] - elev[n] > maxdiff))
                    {
                        maxdiff = elev[c] - elev[n];
                        dmax = d;
                    }
                }*/
                nz = elev[n];
                if ( !processed[n] ) // Compute the spill elevation of the neighbour
                {
                    if( preserve )
					{
						if( nz < (z + mindiff[d]) )
							nz = z + mindiff[d];
					}
					else if( nz <= z )
					{
						nz = z;
                        flowdir[n] = ldd[(d+4)%8];
					}
                    elev[n] = nz;

                    /* if (!flowdir[n] && (elev[c] > elev[n]))
                    {
                        flowdir[n] = d8[(d+4)%8];
                    }
                    elev[n] = std::max(elev[n], elev[c]);*/

                    queue.push(std::move(node(nz, nx, ny)));
                    queued[n] = true;
                }
                /*else if (!flowdir[c]) // Update the steepest gradient direction if needed
                {
                    grad = (z - nz) / length[d];
                    if (grad > maxgrad)
                    {
                        maxgrad = grad;
                        dmax = d;
                    }
                }*/
                
            }
        }
        if (!flowdir[c]) // Record the steepest gradient direction if needed
        {
            flowdir[c] = ldd[getFlowDir(current.x, current.y, z)];
        }
    }

    flowBand->RasterIO(GF_Write, 0, 0, xSize, ySize, flowdir.data(), xSize, ySize, flowBand->GetRasterDataType(), 0, 0);
    spillBand->RasterIO(GF_Write, 0, 0, xSize, ySize, elev, xSize, ySize, spillBand->GetRasterDataType(), 0, 0);

    GDALClose(flowDataset);
    GDALClose(srcDataset);
    GDALClose(spillDataset);
    CPLFree(elev);
    exit(EXIT_SUCCESS);
}