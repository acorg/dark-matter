#!/usr/bin/env sh

pip install -r requirements.txt

# MySQL connector.
tmp=/tmp/dm-install-$$
mkdir $tmp
cd $tmp
curl -O -L http://dev.mysql.com/get/Downloads/Connector-Python/mysql-connector-python-2.1.3.tar.gz
tar xzf mysql-connector-python-2.1.3.tar.gz
cd mysql-connector-python-2.1.3
python setup.py install
