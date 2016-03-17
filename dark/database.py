try:
    from mysql import connector
except ImportError:
    import platform
    if platform.python_implementation() == 'PyPy':
        # PyPy doesn't have a version of mysql connector. Make a fake
        # connector class that raises if used. This allows us to use other
        # 'dark' code that happens to import dark.database but not use it.
        class connector(object):
            @staticmethod
            def connect(**kwargs):
                raise NotImplementedError(
                    'The mysql connector class is not available in PyPy')
    else:
        raise

from os import environ


def getDatabaseConnection():
    return connector.connect(
        host='localhost', user=environ.get('DBI_USER', environ['USER']),
        password=environ['DBI_PASSWORD'], database='ncbi_taxonomy',
        buffered=True)
