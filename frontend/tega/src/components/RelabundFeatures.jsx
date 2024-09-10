import React, { useEffect, useState } from 'react';

const RelabundFeatures = () => {
    const [parameters, setParameters] = useState({
        glomby: '',
        samplesToKeep: '',
        featuresToKeep: '',
        aggregatefeatures: false,
        aggregatefeatures_label: '',
        subsetby: '',
        compareby: '',
        paired: false,
        compareby_order: '',
        colourby: '',
        shapeby: '',
        fillby: '',
        connectby: '',
        facetby: '',
        wrap_facet: false,
        overlay_boxplot: false,
        applyfilters: 'none',
        featcutoff: '',
        GenomeCompletenessCutoff: '',
        PctFromCtgscutoff: '',
        ntop: '',
        minabscorrcoeff: '',
        adjustpval: true,
        padjmeth: 'fdr',
        showonlypbelow: '',
        showonlypadjusted: false,
        maxl2fc: '',
        minl2fc: '',
        addtit: '',
        PPM_normalize_to_bases_sequenced: false,
        log2tran_main_plot: false,
        log2tran_strat_plot: false,
        statsonlog: false,
        y_axis_range: '',
        cdict: '',
        stratify_by_taxlevel: '',
        annotate_phylum: true,
        maxnumplots: '',
        signiflabel: '',
        max_pairwise_cats: 4,
        dump_interpro_descriptions_to_plot: false,
        numthreads: 1,
        nperm: 99,
        ignoreunclassified: true,
        class_to_ignore: 'N_A',
        maxnumtaxa: 20,
        horizontal: true,
        plot_points_on_taxonomy: false,
        use_heatmap_for_stratification: true,
        return_taxon_stratification_df: false,
        return_plots: false,
        rescale_axis_quantiles: '',
    });

    const [relabundfeatureData, setRelabundFeatureData] = useState(null);
    const [objects, setObjects] = useState([]);
    const [filePath, setFilePath] = useState('');
    const [selectedObj, setSelectedObj] = useState('');

    const handleChange = (e) => {
        const { name, value, type, checked } = e.target;
        setParameters({
            ...parameters,
            [name]: type === 'checkbox' ? checked : value
        });
    };

    const handleSubmit = async (e) => {
        e.preventDefault();
        console.log('Selected ExpObj:', selectedObj);
        try {
            // Combine the parameters with selected objects and file path
            const params = {
                filePath,
                ExpObj: selectedObj,
                ...parameters
            };

            // Call IPC Method to run RelabundFeatures script
            const result = await window.electron.runRelabundFeaturesScript(params)

            // Update the AD data with result
            setRelabundFeatureData(result);
        } catch (error) {
            console.error("Error generating AlphaDiversity plot:", error);
        }
    };

    const handleFileUpload = async () => {
        try {
            const filePath = await window.electron.invoke('open-file-dialog');
            if (filePath) {
                setFilePath(filePath);
                const result = await window.electron.invoke('load-rdata-file', filePath);
                setObjects(result);

                // Set default selection to the first object, otherwise it will not default to any
                if (result.length > 0) {
                    setSelectedObj(result[0]);
                }
            }
        } catch (error) {
            console.error('Error opening file dialog:', error);
        }
    };

    const handleObjSelect = (e) => {
        const value = e.target.value;
        setSelectedObj(value);
        console.log('Selected object state:', value);
    };

    const handleClick = () => {
        window.electron.send('navigate-to', 'home');
    };

    const handleDownloadClick = () => {
        // Send IPC event to open the alphadiversity PDF
        window.electron.send('open-alphadiversity-location');
    };

    useEffect(() => {
        // Set default selection if objects array changes
        if (objects.length > 0 && !selectedObj) {
            setSelectedObj(objects[0]);
        }
    }, [objects, selectedObj]);

    return (
        <div>
            <div style={{ position: 'absolute', top: '10px', right: '10px'}}>
                <button onClick={handleClick}>
                    Go Back to Home Page
                </button>
            </div>
        </div>
    );
};

export default RelabundFeatures
