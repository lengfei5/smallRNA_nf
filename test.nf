#!/usr/bin/env nextflow

ch1 = Channel.from(['sample1', 1], ['sample2', 2], ['sample3', 3])
ch2 = Channel.from(['sample1', 4], ['sample3', 6], ['sample2', 5])
ch3 = Channel.from(['sample1', 'A'], ['sample2', 'B'], ['sample3', 'C'])
ch4 = Channel.from(['sample1', 'Z'], ['sample2', 'Y'], ['sample3', 'X'])

ch1.concat(ch2, ch3, ch4)
    .groupTuple()
    .map{ stat1, stat2 -> [stat1, stat2[0], stat2[1], stat2[2], stat2[3]] }
    .println()
