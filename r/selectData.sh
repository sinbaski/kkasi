#!/bin/bash

stocks=(
    AAPL
    ACN
    ADBE
    ADI
    ADP
    ADS
    ADSK
    AKAM
    AMAT
    AVGO
    CA
    CRM
    CSCO
    CTSH
    CTXS
    EA
    EBAY
    EQIX
    FB
    FFIV
    FIS
    FISV
    FSLR
    GOOG
    GOOGL
    HPQ
    HRS
    IBM
    INTC
    INTU
    JNPR
    KLAC
    LLTC
    LRCX
    MA
    MCHP
    MSFT
    MSI
    MU
    NFLX
    NTAP
    NVDA
    ORCL
    PAYX
    PYPL
    QCOM
    QRVO
    RHT
    STX
    SWKS
    SYMC
    TDC
    TEL
    TSS
    TXN
    V
    VRSN
    WDC
    WU
    XLNX
    XRX
    YHOO
);

for i in ${stocks[*]}; do
    # t1=`echo "select min(day), max(day) from ${i}_US;" | mysql -u root -pq1w2e3r4 -N avanza | cut -f 1`;
    # if [[ "$t1" < "2000-01-01" ]]; then
    # 	echo $i
    # fi
    echo "select min(day), max(day) from ${i}_US;"
done | mysql -u root -pq1w2e3r4 -N avanza > /tmp/t1.txt
t1=(`cut -f 1 /tmp/t1.txt`);
t2=(`cut -f 2 /tmp/t1.txt`);

for (( i=1; i <= ${#t1[@]}; i++ )); do
    if [[ "${t1[i]}" < "2000-01-01" && "${t2[i]}" > "2015-01-01" ]]; then
	echo ${stocks[i]}
    fi
done


    
