#!/bin/bash

echo "Get data: "
wget -c http://www.ytzhang.net/files/gait-sp/cache/gaitsp_cache.tar.gz

echo "Extract data:"
tar -xzf gaitsp_cache.tar.gz

echo "Done"

