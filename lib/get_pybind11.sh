#!/bin/sh
  
git clone https://github.com/pybind/pybind11.git pybind11

rm -rf pybind11/{docs,tests,.git,ISSUE_TEMPLATE.md,.gitmodules,.gitignore,.appveyor.yml,.readthedocs.yml}
