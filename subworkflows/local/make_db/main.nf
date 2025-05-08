workflow MAKE_DB {
    take:
    ch_genomes

    main:

    // Update AMRFinderPlus database
    AMRFINDERPLUS_UPDATE()

    // Run MASH sketching
    mash = MASH_SKETCH(ch_genomes)

    // Run AMRFinderPlus with updated DB
    amr = AMRFINDERPLUS_RUN(ch_genomes, AMRFINDERPLUS_UPDATE.out.db)

    // Run MOB-suite
    mob = MOBSUITE_RECON(ch_genomes)

    // Run IntegronFinder
    integrons = INTEGRONFINDER(ch_genomes)

    emit:
    mash_sketches      = mash.out.mash
    amr_reports        = amr.out.report
    mobtyper_results   = mob.out.mobtyper_results
    integron_summaries = integrons.out.summary
}
