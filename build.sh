#!/bin/bash
set -o errexit
PREGEX="^Package: "
PKG=$(grep $PREGEX DESCRIPTION|sed "s/$PREGEX//")
echo Package from DESCRIPTION: $PKG
R -e "if(require(Rcpp))compileAttributes('.')"
R -e "if(require(inlinedocs))package.skeleton.dx('.')"
rm -f man/aum_sort_interface.Rd man/stl_set_insert.Rd man/stl_sort.Rd man/aum_map_interface.Rd  man/aumLineSearch.Rd man/plot.aum_linear_model.Rd man/plot.aum_linear_model_cv.Rd man/predict.aum_linear_model_cv.Rd 
cd ..

RELEASE=$PKG-release
echo Copying $PKG to $RELEASE
rm -rf $RELEASE
cp -r $PKG $RELEASE

echo Editing $RELEASE for CRAN submission
grep -v Remotes $PKG/DESCRIPTION > $RELEASE/DESCRIPTION
rm $RELEASE/tests/testthat/*
cp $PKG/tests/testthat/test-CRAN*.R $RELEASE/tests/testthat

echo Building $RELEASE
RCMD="R --vanilla CMD"
$RCMD build $RELEASE | tee build.out
PKG_TGZ=$(grep building build.out|sed "s/.*\($PKG.*.tar.gz\).*/\1/")

echo Installing $PKG_TGZ
$RCMD INSTALL $PKG_TGZ

echo Checking $PKG_TGZ
$RCMD check --as-cran $PKG_TGZ

echo Checking without any Suggests
R -e "if('check_without_suggests' %in% ls())check_without_suggests('$PKG_TGZ')"
