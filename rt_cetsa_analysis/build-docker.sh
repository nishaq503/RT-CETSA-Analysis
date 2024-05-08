#!/bin/bash

version=$(<VERSION)
docker build . -t polusai/rt_cetsa_analysis_test_7:${version}