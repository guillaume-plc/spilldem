#include <getopt.h>
#include <iostream>
#include <queue>
#include "gdal_priv.h"
#include "cpl_conv.h"

#include "SpillDEM.h" // config file

static void usage(const char* name)
{
    printf("%s version %d.%d\n"
           "usage: %s <options> datasource\n"
           "Options:\n"
            "\t-o, --output        output file name\n"
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
        {"verbose", no_argument, nullptr, 'v'},
        {"help", no_argument, nullptr, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    bool verbose = false;
    std::string in_file = "";
    std::string out_file = "out.tif";
    while ((opt = getopt_long(argc, argv, ":o:vh", long_opts, nullptr)) != -1) 
    {
        switch (opt) 
        {
        case 'v':
            verbose = true;
            break;
        case 'o':
            out_file = std::string(optarg);
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
    in_file = argv[optind];

    GDALAllRegister();
    const char *format = "GTiff";
    GDALDriver *driver;
    driver = GetGDALDriverManager()->GetDriverByName(format);
    if( driver == nullptr )
    {
        exit(EXIT_FAILURE);
    }

    GDALDataset *srcDataset;
    srcDataset = (GDALDataset *)GDALOpen(in_file.c_str(), GA_ReadOnly);
    if ( srcDataset == nullptr)
    {
        exit(EXIT_FAILURE);
    }

    GDALDataset *dstDataset;
    dstDataset = driver->CreateCopy(out_file.c_str(), srcDataset, FALSE,
                                      NULL, NULL, NULL );
    if ( dstDataset == nullptr )
    {
        GDALClose(srcDataset);
        exit(EXIT_FAILURE);
    }

    GDALRasterBand *srcBand, *dstBand;
    srcBand = srcDataset->GetRasterBand(1);
    dstBand = dstDataset->GetRasterBand(1);
    const int nx = srcBand->GetXSize(), ny = srcBand->GetYSize();
    double nodata = srcBand->GetNoDataValue();

    float *dem, *filleddem;
    dem = (float *) CPLMalloc(sizeof(float)*nx*ny);
    srcBand->RasterIO(GF_Read, 0, 0, nx, ny, dem, nx, ny, GDT_Float32, 0, 0);
    dstBand->Fill(nodata);

    filleddem = (float *) CPLMalloc(sizeof(float)*nx*ny);
    std::array<dir, 8> ngh = {
        dir(1, 0),
        dir(1, -1),
        dir(0, -1),
        dir(-1, -1),
        dir(-1, 0),
        dir(-1, 1),
        dir(0, 1),
        dir(1, 1)
    };

    std::array<unsigned char, 8> d8 = {1, 2, 4, 8, 16, 32, 64, 128};
    std::priority_queue<node> queue;
    std::vector<bool> queued(nx*ny, false);
    std::vector<bool> processed(nx*ny, false);
    for (int x = 0; x < nx; x++) // fill edge cell with the elevation value
    {
        for (int y = 0; y < ny; y++)
        {
            bool edge;
            if (dem[y * nx + x] == nodata)
            {
                edge = false;
                processed[y * nx + x] = true;
            }
            else if (x == 0 || x == nx-1 || y == 0 || y == ny-1)
            {
                edge = true;
            }
            else
            {
                edge = false;
                for (int d = 0; d < 8; d++)
                {
                    if (dem[(y + ngh[d].dy) * nx + (x + ngh[d].dx)] == nodata)
                    {
                        edge = true;
                        break;
                    }
                }
            }
            filleddem[y * nx + x] = edge ? dem[y * nx + x] : nodata;
            if (edge)
            {
                queue.push(node(filleddem[y * nx + x], x, y));
                queued[y * nx + x] = true;
            }
        }
    }

    while (!queue.empty())
    {
        node c(queue.top());
        queue.pop();
        int cid = c.y*nx + c.x;
        processed[cid] = true;
        queued[cid] = false;
        for (int d = 0; d < 8; d++)
        {
            int nghy = c.y+ngh[d].dy, nghx = c.x+ngh[d].dx; 
            int nid = nghy * nx + nghx;
            if (nghy < 0 || nghy >= ny || nghx < 0 || nghx >= nx || queued[nid] || processed[nid])
            {
                continue;
            }
            else
            {
                filleddem[nid] = std::max(dem[nid], filleddem[cid]);
                queue.push(node(filleddem[nid], nghx, nghy));
                queued[nid] = true;
            }
        }
    }

    dstBand->RasterIO(GF_Write, 0, 0, nx, ny, filleddem, nx, ny, dstBand->GetRasterDataType(), 0, 0);

    CPLFree(dem);
    CPLFree(filleddem);
    GDALClose(dstDataset);
    GDALClose(srcDataset);
    exit(EXIT_SUCCESS);
}