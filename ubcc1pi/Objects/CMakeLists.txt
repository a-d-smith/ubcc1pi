art_make(
    LIBRARY_NAME
        ubcc1pi_Objects
    LIB_LIBRARIES
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_PRINCIPAL}
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES_BASIC}
        ${ART_PERSISTENCY_COMMON}
        ${ART_PERSISTENCY_PROVENANCE}
        ${ROOT_BASIC_LIB_LIST}
        canvas
        cetlib cetlib_except
        nusimdata_SimulationBase
        lardataobj_RecoBase
        lardataobj_AnalysisBase
        art_Framework_Principal
        art_Persistency_Common
        ubcc1pi_Helpers
        ubcc1pi_Interface
    )

install_headers()
install_source()
install_fhicl()
