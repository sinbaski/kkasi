#!/bin/bash

symbols=`echo "select symbol from SP500_components;" | mysql -N -u"sinbaski" -p"q1w2e3r4" -h localhost avanza`;
# symbols=`echo "select symbol from OMXS30_components;" | mysql -N -u"sinbaski" -p"q1w2e3r4" -h localhost avanza`;

for s in `echo $symbols`; do
    t=`echo $s | sed -e s/[.]/_/g | sed -e s/[-]/_series_/g`;
    echo "CREATE TABLE $t-US ("
    echo "  day date NOT NULL,"
    echo "  high decimal(7,2) NOT NULL,"
    echo "  low decimal(7,2) NOT NULL,"
    echo "  closing decimal(7,2) NOT NULL,"
    echo "  volume int(11) DEFAULT NULL,"
    echo "  PRIMARY KEY (day)"
    echo ");"
done

# | mysql -u sinbaski -pq1w2e3r4 avanza

