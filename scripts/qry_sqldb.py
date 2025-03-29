#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
qry_sqldb.py

This script connects to a PostgreSQL database, executes a SQL query
to retrieve  data in the ChEMBL database, and saves the resulting
data to a CSV file.

Input:
- SQL query as a text file

Output:
- CSV file

Dependencies:
- sys
- argparse
- pandas
- adbc_driver_postgresql
- dotenv

Usage:
- Ensure that the .env file contains the correct database URI
  under the key 'DB_URI'
- Run the script using Python 3.

Author: Brad Dallin
Date: March 4, 2025
"""


########################################################################
## Imports
########################################################################
import sys
import argparse
import pandas as pd
import adbc_driver_postgresql.dbapi
from dotenv import dotenv_values


########################################################################
## Functions
########################################################################
def parse_args(argv) -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        prog="qry_sqldb.py",
        description="Retreive data from ChEMBL PostgreSQL and save to a CSV file",
    )
    parser.add_argument(
        "input",
        action="store",
        type=str,
        help="Text file containing SQL query",
    )
    parser.add_argument(
        "output",
        action="store",
        type=str,
        help="CSV file with queried data",
    )
    parser.add_argument(
        "--env_path",
        action="store",
        type=str,
        default=".env",
        required=False,
        help="Path to .env file"
    )
    args = parser.parse_args(argv)
    return args


# Read query
def read_query(file_path: str) -> str:
    """Read and store SQL query from text file"""
    with open(file_path, "r",
              encoding="utf-8") as txt_file:
        txt = txt_file.read()
    return txt


# Connect to PostgreSQL database
def sql_query(
        qry: str,
        env_path: str = ".env"
        ) -> pd.DataFrame:
    """Run SQL query and return as Pandas DataFrame"""
    config = dotenv_values(env_path)
    # Test connection
    try:
        uri = config["DB_URI"]
        conn = adbc_driver_postgresql.dbapi.connect(uri)
        with conn.cursor() as cur:
            cur.execute("SELECT 1")
            assert cur.fetchone() == (1,)
        print("Database connected successfully\n")
    except adbc_driver_postgresql.dbapi.DatabaseError as err:
        print(f"Database not connected successfully: {err}\n")
        sys.exit(1)

    # Run input qry
    with conn.cursor() as cur:
        cur.execute(qry)
        df = pd.DataFrame(
            cur.fetchall(),
            columns=[desc[0] for desc in cur.description]
        )
    print(f"{df.shape[0]} molecules found!\n")
    return df


# Main function
def main(argv):
    """Main function"""
    args = parse_args(argv)
    sql = read_query(args.input)
    df = sql_query(sql, args.env_path)
    # Save to file
    df.to_csv(
        args.output,
        index=False,
        encoding="utf-8"
    )
    return 0


########################################################################
## Run
########################################################################
if __name__ == "__main__":
    main(sys.argv[1:])
