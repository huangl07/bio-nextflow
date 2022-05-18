#!/usr/bin/env nextflow
params.worksh = 'fq.list'
params.help = false
params.cpu = 1 
params.mem = "3G"
params.mode = "split"
params.excutor = "local"
params.queue = "DNA"
params.outdir="./"
def helpMessage() {
   
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --worksh worksh

    --worksh <file> work sh
    --cpu   <num>   threads cpu number
    --mem   <str>   mem resource  default 3G
    --mode  default split 
            submit for submit a shell for a single node 
            split  for qsub sum shell for parallel 
    --excutor   default slurm 
            slurm
            sge
            local
    --queue     queue default DNA for 10.2.4.236
    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}

if (params.mode == "split"){
    jobs = file(params.worksh).readLines()
    def count = 1..jobs.size()
    process qsub{
        publishDir "log" , pattern: "*log"
        queue params.queue
        cpus params.cpu
        memory params.mem
        executor params.excutor
        input:
            val(sh) from jobs
            val(num) from count
        output:
            file "*.log"
        script:
        """
            $sh 2>${num}.log 1>${num}.log
            echo "$sh" >>${num}.log
        """
    }
}else{
    process submit{
        publishDir "log" , pattern: "*log"
        queue params.queue
        cpus params.cpu
        memory params.mem
        executor params.excutor
        input:
            file sh from file(params.worksh)
        output:
            file  "run.log"
        """
            sh ${sh} 2>run.log 1>run.log
        """
    }
}
workflow.onComplete {
    println "Pipeline completed!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

