#!/bin/bash

protoc --cpp_out=. ./vector_2d.proto
protoc --cpp_out=. ./range_2d.proto
protoc --cpp_out=. ./grid_2d.proto
protoc --cpp_out=. ./face_array_2d.proto
protoc --cpp_out=. ./face_array_ap_2d.proto
