#!/bin/sh

set -e
set -x

cp plugins/MultiScaleTubularityMeasure_Plugin.jar ~/fiji/plugins/
rm -f ~/fiji/plugins/ITK/linux/*
cp ../../ITK/bin/*.so ~/fiji/plugins/ITK/linux/
cp build/linux/libMultiScaleTubularityMeasure.so ~/fiji/plugins/ITK/linux/
