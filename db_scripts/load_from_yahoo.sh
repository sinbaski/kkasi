#!/bin/bash

symbols=`echo "select symbol from SP500_components where available=1;" | mysql -N -u"sinbaski" -p"q1w2e3r4" -h localhost avanza`;
# symbols=`echo "select symbol from OMXS30_components;" | mysql -N -u"sinbaski" -p"q1w2e3r4" -h localhost avanza`;
# symbols=`echo "select symbol from indices;" | mysql -N -u"sinbaski" -p"q1w2e3r4" -h localhost avanza`;

for symbol in `echo $symbols`; do
    # if test -f ~/Downloads/$symbol.csv; then
    # 	rm ~/Downloads/$symbol.csv
    # fi
    if [ -z "$str" ]; then str=$symbol; else str="$str,$symbol"; fi
done
curl http://real-chart.finance.yahoo.com/table.csv?s={"$str"} -o /tmp/#1.csv
if [ ! -d /tmp/sql ]; then
    mkdir /tmp/sql
fi
# python ./load_from_yahoo.py
    # tbl=`echo $symbol | sed -e 's:\.:_:g'`;
    # tbl=`echo $symbol | sed -e 's:\.:_:g' | sed -r -e 's:$:_US:g'`;
#     echo "drop table if exists $tbl;" > /tmp/$symbol.txt;
#     echo "create table if not exists $tbl (
# day date primary key,
# high decimal(7,2) not null,
# low decimal(7,2) not null,
# closing decimal(7,2) not null,
# volume int
# );" >> /tmp/$symbol.txt;

#     if [ -f /tmp/$symbol.txt ]; then rm /tmp/$symbol.txt; fi
#     echo "insert into $tbl values " >> /tmp/$symbol.txt;
#     tail -n +2 ~/Downloads/$symbol.csv |
#     head -n 7 |
#     cut -f 1,3,4,5,6 -d',' |
#     awk -F, '
# {
# if (FNR == 1) printf "(\"%s\", %.2f, %.2f, %.2f, %d)", $1,$2,$3,$4,$5
# else printf ",\n(\"%s\", %.2f, %.2f, %.2f, %d)", $1,$2,$3,$4,$5
# }
# END {print ";"}' >> /tmp/$symbol.txt;


    # printf "drop table if exists %s;" $symbol > /tmp/$symbol.txt
    # mysql -h 83.176.196.41 -u sinbaski -p'q1w2e3r4' avanza < /tmp/$symbol.txt;

    # mysql -h localhost -u sinbaski -p'q1w2e3r4' avanza < /tmp/$symbol.txt;


