from os import environ

from mysql import connector


def getDatabaseConnection():
    return connector.connect(
        host='localhost', user=environ.get('DBI_USER', environ['USER']),
        password=environ['DBI_PASSWORD'], database='ncbi_taxonomy',
        buffered=True)
