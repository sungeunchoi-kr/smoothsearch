#!/bin/bash

while true
do
    START=$(date +%s%3N);
    block_addr=$(curl -s http://162.243.237.110:8718/v1/free-blocks)
    block_addr="${block_addr%\"}"
    block_addr="${block_addr#\"}"
    echo "got block address $block_addr."

    a=$(./run $block_addr)
    echo "$a"

    END=$(date +%s%3N);
    echo "65536 block searched in $((END-START)) ms."

    curl -s -d "$a" "http://162.243.237.110:8718/v1/searched-blocks/$block_addr"
    echo ""
done

