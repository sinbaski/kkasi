#!/usr/bin/python

import os
import sys
import _mysql
import csv

# conn = _mysql.connect(host=os.environ.get("PB"),
#                      user="sinbaski",
#                      passwd="q1w2e3r4",
#                      db="avanza");
conn = _mysql.connect(host="localhost",
                     user="sinbaski",
                     passwd="q1w2e3r4",
                     db="avanza");
index = "SP500";
suffix = "_US";
conn.query("select symbol from {0}_components where available=1;".format(index));


Symbols = conn.store_result();
while True:
    symbol = Symbols.fetch_row(maxrows=1, how=0);
    if not symbol:
        break;

    tbl = symbol[0][0].replace(".", "_");
    tbl = tbl.replace("-", "_series_");
    tbl += suffix;
    stmt = 'insert into {0} values'.format(tbl)
    csv_file_name = '/tmp/{0}.csv'.format(symbol[0][0]);
    
    conn.query('select max(day) as t from {0};'.format(tbl));
    X = conn.use_result().fetch_row(maxrows = 1);
    # last_date = X[0][0];
    # conn.query('delete from {0};'.format(tbl));

    num = 0;
    with open(csv_file_name, 'rb') as csvfile:
        my_reader = csv.reader(csvfile)
        c = 0;
        for row in my_reader:
            c += 1;
            # Skip the header
            if c == 1 or row[2] == row[3]:
                continue
            elif X != None and row[0] <= X[0][0]:
                break;
            num += 1;
            if num == 1:
                stmt += ' ("{0}", {1}, {2}, {3}, {4})'.format(
                    row[0], row[2], row[3], row[4], row[5]);
            else:
                stmt += ', ("{0}", {1}, {2}, {3}, {4})'.format(
                    row[0], row[2], row[3], row[4], row[5]);
    print '%d rows to insert into %s.' % (num, tbl);
    if (num):
        # stmt += ' on duplicate key update ';
        # stmt += 'high=values(high), low=values(low),';
        # stmt += 'closing=values(closing), volume=values(volume)';
        stmt += ';';
        logfile = open('/tmp/sql/%s.sql' % symbol[0][0], 'wb');
        logfile.write('%s\n' % stmt);
        logfile.close();
        conn.query(stmt);
    
conn.close();

