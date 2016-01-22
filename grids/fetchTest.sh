#!/bin/bash
# Fetch APPLgrids for testing the APFELgrid plugin
wget http://www.hepforge.org/archive/applgrid/atlas-Z0-arxiv-1109.5141.tgz
tar -xvzf ./atlas-Z0-arxiv-1109.5141.tgz --strip-components 2
rm atlas-Z0-arxiv-1109.5141.tgz