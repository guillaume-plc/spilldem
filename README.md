# spillDEM
Fast DEM surface depressions filling using the spill elevation and least-cost search approach from [Wang, Lei & Liu, Holiday. (2006)](http://dx.doi.org/10.1080/13658810500433453).

Basically a standalone version of SAGA's implementation, relying only on GDAL for IO (see [SAGA's source code](https://sourceforge.net/projects/saga-gis/) for the 'Fill sinks (Wang & Liu)' preprocessing algorithm). As for SAGA's implementation the original algorithm is modified to allow for preservation of a minimum slope gradient between cells.

## Requirements
- GDAL v3+

## Installation

## Usage
