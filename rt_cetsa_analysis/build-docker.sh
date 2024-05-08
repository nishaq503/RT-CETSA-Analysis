#!/bin/bash

version=$(<VERSION)
docker build . -t polusai/rt_cetsa_analysis:${version}