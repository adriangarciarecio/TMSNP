#!/usr/bin/env python3

################################################
import sqlalchemy as sql
import pandas as pd

# Set the connection
# engine = sql.create_engine(
#     f"mysql://lmcdb:{os.getenv('LMCDB_PASS')}@alf03.uab.cat/tmsnp"
# )
engine = sql.create_engine(f"mysql://adrian:compmod5@localhost/tmsnp")
connection = engine.connect()

connection.execute("DROP TABLE IF EXISTS acc_code")
connection.execute(
    """CREATE TABLE acc_code (acc varchar(20), id varchar(30), gene varchar(20));"""
)

file_data = open("./ensembl_list.txt", "r")
for line in file_data:
    l = line.split("\t")
    connection.execute(
        "INSERT INTO acc_code VALUES (%s, %s, %s)", (l[0], l[1], l[2][:-1])
    )
    print(f"> Updating the info of {l[0]}.")
    connection.execute(f'UPDATE snps SET id = "{l[1]}" WHERE snps.acc LIKE "{l[0]}";')
    connection.execute(
        f'UPDATE snps SET gene = "{l[2][:-1]}" WHERE snps.acc LIKE "{l[0]}";'
    )

connection.close()

