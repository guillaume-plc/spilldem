#include <getopt.h>
#include <iostream>
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
    const char *pszFormat = "GTiff";
    GDALDriver *poDriver;
    poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    if( poDriver == nullptr )
    {
        exit(EXIT_FAILURE);
    }

    GDALDataset *dem;
    dem = (GDALDataset *)GDALOpen(in_file.c_str(), GA_ReadOnly);
    if ( dem == nullptr)
    {
        exit(EXIT_FAILURE);
    }
    GDALDataset *filled_dem;
    filled_dem = poDriver->CreateCopy(out_file.c_str(), dem, FALSE,
                                      NULL, NULL, NULL );
    if ( filled_dem == nullptr)
    {
        GDALClose(dem);
        exit(EXIT_FAILURE);
    }

    GDALClose(filled_dem);
    GDALClose(dem);
    exit(EXIT_SUCCESS);
}