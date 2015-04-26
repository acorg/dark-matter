import mysql.connector

from os import environ


def getDatabaseConnection():
    return mysql.connector.connect(
        host='localhost', user=environ.get('DBI_USER', environ['USER']),
        password=environ['DBI_PASSWORD'], database='ncbi_taxonomy')
