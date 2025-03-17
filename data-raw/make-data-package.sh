#!/bin/bash

PROJECT_DIR=~/ebirdst/ebird-best-practices/
PKG_DIR=${PROJECT_DIR}/ebird-best-practices-data/
rm -rf $PKG_DIR
mkdir -p ${PKG_DIR}/data/
mkdir -p ${PKG_DIR}/data-raw/

# cp files into package folder
cp ${PROJECT_DIR}/data/checklists-zf_woothr_jun_us-ga.csv ${PKG_DIR}/data/
cp ${PROJECT_DIR}/data/environmental-variables_checklists_jun_us-ga.csv ${PKG_DIR}/data/
cp ${PROJECT_DIR}/data/environmental-variables_prediction-grid_us-ga.csv ${PKG_DIR}/data/
cp ${PROJECT_DIR}/data/gis-data.gpkg ${PKG_DIR}/data/
cp ${PROJECT_DIR}/data/prediction-grid_us-ga.tif ${PKG_DIR}/data/
cp ${PROJECT_DIR}/data-raw/elevation_gmted_1km_us-ga.tif ${PKG_DIR}/data-raw/
cp ${PROJECT_DIR}/data-raw/landcover_mcd12q1_umd_us-ga_2014-2022.tif ${PKG_DIR}/data-raw/
cp ${PROJECT_DIR}/data-raw/mcd12q1_umd_classes.csv ${PKG_DIR}/data-raw/

# compress
cd $PROJECT_DIR
zip -rv ebird-best-practices-data.zip ebird-best-practices-data/  -x "*/.*" ".*"
mv ebird-best-practices-data.zip ${PROJECT_DIR}/data-raw/
cd -

rm -rf $PKG_DIR
