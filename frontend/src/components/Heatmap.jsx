import React, { useEffect, useState } from 'react';
import Button from '@mui/material/Button';
import CloudUploadIcon from '@mui/icons-material/CloudUpload';

const Heatmap = ({ handleNavigateTo }) => {
    const [parameters, setParameters] = useState({
        glomby: '',
        hmtype: 'exploratory',
        samplesToKeep: '',
        featuresToKeep: '',
        applyfilters: 'none',
        featcutoff: '',
        GenomeCompletenessCutoff: '',
        PPM_normalize_to_bases_sequenced: false,
        subsetby: '',
        compareby: '',
        invertbinaryorder: false,
        showonlypbelow: '',
        adj_pval_for_threshold: false,
        ntop: '',
        colcategories: '',
        splitcolsby: '',
        ordercolsby: '',
        textby: '',
        cluster_samples_per_heatmap: false,
        cluster_features_per_heatmap: false,
        label_samples: false,
        cluster_rows: false,
        max_rows_in_heatmap: '',
        no_underscores: false,
        showGram: false,
        show_GenomeCompleteness: false,
        addtit: '',
        hmasPA: false,
        threshPA: 0,
        cluster_column_slices: true,
        column_split_group_order: '',
        row_order: '',
        discard_SDoverMean_below: '',
        maxl2fc: '',
        minl2fc: '',
        fun_for_l2fc: 'geom_mean',
        showpval: true,
        showroundedpval: true,
        showl2fc: true,
        assay_for_matrix: 'BaseCounts',
        normalization: 'relabund',
        asPPM: true,
        scaled: false,
        cdict: '',
        maxnumheatmaps: '',
        numthreads: 1,
        statsonlog: false,
        ignoreunclassified: true,
        returnstats: false,
        class_to_ignore: 'N_A'
    });

    const displayNames = {
        glomby: 'Glom By',
        hmtype: 'Heatmap Type',
        samplesToKeep: 'Samples to Keep',
        featuresToKeep: 'Features to Keep',
        applyfilters: 'Apply Filters',
        featcutoff: 'Feature Cutoff',
        GenomeCompletenessCutoff: 'Genome Completeness Cutoff',
        PPM_normalize_to_bases_sequenced: 'PPM Normalize to Bases Sequenced',
        subsetby: 'Subset By',
        compareby: 'Compare By',
        invertbinaryorder: 'Invert Binary Order',
        showonlypbelow: 'Show Only P Below',
        adj_pval_for_threshold: 'Adjust P-Value for Threshold',
        ntop: 'N Top',
        colcategories: 'Column Categories',
        splitcolsby: 'Split Columns By',
        ordercolsby: 'Order Columns By',
        textby: 'Text By',
        cluster_samples_per_heatmap: 'Cluster Samples per Heatmap',
        cluster_features_per_heatmap: 'Cluster Features per Heatmap',
        label_samples: 'Label Samples',
        cluster_rows: 'Cluster Rows',
        max_rows_in_heatmap: 'Max Rows in Heatmap',
        no_underscores: 'No Underscores',
        showGram: 'Show Gram',
        show_GenomeCompleteness: 'Show Genome Completeness',
        addtit: 'Add Title',
        hmasPA: 'Heatmap as PA',
        threshPA: 'Threshold PA',
        cluster_column_slices: 'Cluster Column Slices',
        column_split_group_order: 'Column Split Group Order',
        row_order: 'Row Order',
        discard_SDoverMean_below: 'Discard SD over Mean Below',
        maxl2fc: 'Max L2FC',
        minl2fc: 'Min L2FC',
        fun_for_l2fc: 'Function for L2FC',
        showpval: 'Show P-Value',
        showroundedpval: 'Show Rounded P-Value',
        showl2fc: 'Show L2FC',
        assay_for_matrix: 'Assay for Matrix',
        normalization: 'Normalization',
        asPPM: 'As PPM',
        scaled: 'Scaled',
        cdict: 'Color Dictionary',
        maxnumheatmaps: 'Max Number of Heatmaps',
        numthreads: 'Number of Threads',
        statsonlog: 'Stats on Log',
        ignoreunclassified: 'Ignore Unclassified',
        returnstats: 'Return Stats',
        class_to_ignore: 'Class to Ignore'
    };

    const [heatmapData, setHeatmapData] = useState(null);
    const [objects, setObjects] = useState([]); // for ExpObj dropdown
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
                advancedSettings: {},
                ...parameters
            };

            // Call IPC method to run the heatmap script
            const result = await window.electron.runHeatmapScript(params);

            // Update the heatmap data with result
            setHeatmapData(result);
        } catch (error) {
            console.error('Error generating heatmap:', error);
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
