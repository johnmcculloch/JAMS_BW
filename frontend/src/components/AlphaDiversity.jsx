import React, { useEffect, useState } from 'react';
import Button from '@mui/material/Button';
import CloudUploadIcon from '@mui/icons-material/CloudUpload';

const AlphaDiversity =({ handleNavigateTo }) => {
    const [parameters, setParameters] = useState({
        measures: 'c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", "GeneCount")',
        stratify_by_kingdoms: true,
        glomby: '',
        samplesToKeep: '',
        subsetby: '',
        compareby: '',
        compareby_order: '',
        colourby: '',
        shapeby: '',
        fillby: '',
        pairby: '',
        connectby: '',
        facetby: '',
        wrap_facet: '',
        overlay_boxplot: false,
        applyfilters: 'none',
        featcutoff: '',
        GenomeCompletenessCutoff: '',
        PPM_normalize_to_bases_sequenced: false,
        addtit: '',
        signiflabel: 'p.format',
        max_pairwise_cats: '4',
        ignoreunclassified: true,
        class_to_ignore: 'N_A',
    });
    
}
