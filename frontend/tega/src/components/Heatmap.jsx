import React, { useEffect, useState } from 'react';
import Button from '@mui/material/Button';
import CloudUploadIcon from '@mui/icons-material/CloudUpload';


const Heatmap = () => {
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

    const handleObjSelect = (e) => {
        const value = e.target.value;
        setSelectedObj(value);
        console.log('Selected object state:', value);
    };


    const handleClick = () => {
        window.electron.send('navigate-to', 'home');
    };

    const handleDownloadClick = () => {
        // Send IPC event to open the heatmap PDF
        window.electron.send('open-heatmap-location');
    };

    useEffect(() => {
        // Set default selection if objects array changes
        if (objects.length > 0 && !selectedObj) {
            setSelectedObj(objects[0]);
        }
    }, [objects, selectedObj]);

    return (
        <div>
            <div style={{ position: 'absolute', top: '10px', right: '10px' }}>
                <button onClick={handleClick}>
                    Go Back to Home Page
                </button>
            </div>


            <h1>Generate Heatmap</h1>
        <div>

            {/* File upload for RData file */}
            <h3>Upload R Data File for Heatmap</h3>
            <Button 
            component="label"
            variant="contained" 
            startIcon={<CloudUploadIcon />}
                >
                    Upload RData File
                    <input
                    type='file'
                    accept='.rdata, .rda'
                    onChange={handleFileUpload}
                    style={{ display: 'none' }}
                    />
                    </Button>
        </div>

            {/* Dropdown for selecting Summarized Experiment Object */}
            {objects.length > 0 && (
                <div>
                    <h3>Select Summarized Experiment Object</h3>
                    <select onChange={handleObjSelect} value={selectedObj}>
                        {objects.map((obj, index) => (
                            <option key={index} value={obj}>
                                {obj}
                            </option>
                        ))}
                    </select>
                </div>
            )}


            {/* Form for Heatmap Parameters */}
            <form onSubmit={handleSubmit}>
                {Object.keys(parameters).map((key) => (
                    <div key={key}>
                        <label htmlFor={key}>{key}:</label>
                        {key === 'applyfilters' ? (
                            <select
                            id={key}
                            name={key}
                            value={parameters[key]}
                            onChange={handleChange}
                            >
                                <option value='none'>none</option>
                                <option value='light'>light</option>
                                <option value='moderate'>moderate</option>
                                <option value='stringent'>stringent</option>
                            </select>
                        ) : typeof parameters[key] === 'boolean' ? (
                            <input
                                type="checkbox"
                                id={key}
                                name={key}
                                checked={parameters[key]}
                                onChange={handleChange}
                            />
                        ): (
                            <input
                            type="text"
                            id={key}
                            name={key}
                            value={parameters[key]}
                            onChange={handleChange}
                            />
                        )}
                    </div>
                ))}
                <Button type='submit' variant='outlined' color='primary'>
                    Generate Heatmap
                </Button>
            </form> 

            {heatmapData && (
                <div id="heatmap=container">
                    <h2>Heatmap</h2>
                    <Button
                    variant="contained"
                    color='primary'
                    onClick={handleDownloadClick}
                    >
                        Open Heatmap PDF
                    </Button>
                </div>
            )}
        </div>
    );
};

export default Heatmap;
