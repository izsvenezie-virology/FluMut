from peewee import DatabaseProxy

from flumut.flumutdb.initializer import initialize

DATABASE_PROXY = DatabaseProxy()

initialize()
