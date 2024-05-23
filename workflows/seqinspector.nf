/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                        } from '../modules/nf-core/fastqc/main'

include { MULTIQC as MULTIQC_GLOBAL     } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_PER_LANE   } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_PER_GROUP  } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_PER_RUNDIR } from '../modules/nf-core/multiqc/main'

include { paramsSummaryMap              } from 'plugin/nf-validation'
include { paramsSummaryMultiqc          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText        } from '../subworkflows/local/utils_nfcore_seqinspector_pipeline'
include { SORT_FILES as SORT_LANES      } from '../subworkflows/local/utils_nfcore_seqinspector_pipeline'
include { SORT_FILES as SORT_GROUPS     } from '../subworkflows/local/utils_nfcore_seqinspector_pipeline'
include { SORT_FILES as SORT_RUNDIRS    } from '../subworkflows/local/utils_nfcore_seqinspector_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SEQINSPECTOR {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions            = Channel.empty()
    ch_multiqc_files       = Channel.empty()
    ch_multiqc_extra_files = Channel.empty()
    ch_multiqc_reports     = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_logo   = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params                        = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(
        paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_extra_files = ch_multiqc_extra_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_extra_files = ch_multiqc_extra_files.mix(ch_collated_versions)
    ch_multiqc_extra_files = ch_multiqc_extra_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC_GLOBAL (
        ch_multiqc_files
            .map { meta, file -> file }
            .mix(ch_multiqc_extra_files)
            .collect(),
        ch_multiqc_config.toList(),
        Channel.empty().toList(),
        ch_multiqc_logo.toList()
    )

    // Generate reports by lane
    (sorted_samples_by_lane, config_by_lane) = SORT_LANES(
        ch_multiqc_files.toList(),
        ch_multiqc_extra_files.toList(),
        "lane",
    )

    MULTIQC_PER_LANE(
        sorted_samples_by_lane,
        ch_multiqc_config.toList(),
        config_by_lane,
        ch_multiqc_logo.toList()
    )

    // Generate reports by group
    (sorted_samples_by_group, config_by_group) = SORT_GROUPS(
        ch_multiqc_files.toList(),
        ch_multiqc_extra_files.toList(),
        "group",
    )

    MULTIQC_PER_GROUP(
        sorted_samples_by_group,
        ch_multiqc_config.toList(),
        config_by_group,
        ch_multiqc_logo.toList()
    )

    // Generate reports by rundir
    (sorted_samples_by_rundir, config_by_rundir) = SORT_RUNDIRS(
        ch_multiqc_files.toList(),
        ch_multiqc_extra_files.toList(),
        "rundir",
    )

    MULTIQC_PER_RUNDIR(
        sorted_samples_by_rundir,
        ch_multiqc_config.toList(),
        config_by_rundir,
        ch_multiqc_logo.toList()
    )


    emit:
    global_report = MULTIQC_GLOBAL.out.report.toList()             // channel: /path/to/multiqc_report.html
    lane_reports = MULTIQC_PER_LANE.out.report.toList()     // channel: [ /path/to/multiqc_report.html ]
    group_reports = MULTIQC_PER_GROUP.out.report.toList()   // channel: [ /path/to/multiqc_report.html ]
    rundir_reports = MULTIQC_PER_RUNDIR.out.report.toList() // channel: [ /path/to/multiqc_report.html ]
    versions       = ch_versions                            // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
