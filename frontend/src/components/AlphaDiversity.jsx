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

    const displayNames = {
        measures: 'Measures',
        stratify_by_kingdoms: 'Stratify by Kingdoms',
        glomby: 'Glom By',
        samplesToKeep: 'Samples to Keep',
        subsetby: 'Subset By',
        compareby: 'Compare By',
        compareby_order: 'Order to Compare By',
        colourby: 'Colour By',
        shapeby: 'Shape By',
        fillby: 'Fill By',
        pairby: 'Pair By',
        connectby: 'Connect By',
        facetby: 'Facet By',
        wrap_facet: 'Wrap Facet',
        overlay_boxplot: 'Overlay Boxplot',
        applyfilters: 'Apply Filters',
        featcutoff: 'Feature Cutoff',
        GenomeCompletenessCutoff: 'Genome Completeness Cutoff',
        PPM_normalize_to_bases_sequenced: 'PPM Normalize to Bases Sequenced',
        addtit: 'Add Title',
        signiflabel: 'Significance Label',
        max_pairwise_cats: 'Maximum Pairwise Categories',
        ignoreunclassified: 'Ignore Unclassified',
        class_to_ignore: 'Class to Ignore'
    };

    const [alphaDiversityData, setAlphaDiversityData] = useState(null);
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

            // Call IPC method to run AlphaDiversity script
            const result = await window.electron.runAlphaDiversityScript(params)

            // Update the AD data with result
            setAlphaDiversityData(result);
        } catch (error) {
            console.error("Error generating AlphaDiversity plot:", error)
        }
    };

    const handleFileUpload = async () => {
        try {
            const filePath = await window.electron.invoke('open-file-dialog');
            if (filePath) {
                setFilePath(filePath);
                const result = await window.electron.invoke('load-rdata-file', filePath);
                setObjects(result);

                // Set default selection to the first object, otherwise it won't default to any
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


}
