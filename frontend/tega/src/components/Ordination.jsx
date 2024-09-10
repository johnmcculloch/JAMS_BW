import React, { useEffect, useState } from 'react';
import Button from '@mui/material/Button';
import CloudUploadIcon from '@mui/icons-material/CloudUpload';


const Ordination = () => {
    const [parameters, setParameters] = useState({
        glomby: '',
        algorithm: '',
        PCA_Components: '',
        distmehod: '',
        subsetby: '',
        compareby: '',
        colourby: '',
        colorby: '',
        shapeby: '',
        samplesToKeep: '',
        featuresToKeep: '',
        samplesToHighlight: '',
        ignoreunclassified: true,
        applyfilters: 'none',
        featcutoff: '',
        GenomeCompletenessCutoff: '',
        discard_SDoverMean_below: '',
        asPPM: true,
        normalization: 'relabund',
        PPM_normalize_to_bases_sequenced: false,
        assay_for_matrix: 'BaseCounts',
        use_letters_as_shapes: false,
        sizeby: '',
        connectby: '',
        connection_orderby: '',
        textby: '',
        ellipseby: '',
        dotsize: '',
        log2tran: false,
        tsne_perplx: '',
        max_neighbors: '',
        permanova: false,
        plotcentroids: false,
        highlight_centroids: false,
        show_centroid_distances: false,
        calculate_centroid_distances_in_all_dimensions: false,
        addtit: '',
        cdict: '',
        grid: false,
        forceaspectratio: '',
        numthreads: '',
        return_coordinates_matrix: false,
        permanova_permutations: '',
        include_components_variance_plot: false,
        class_to_ignore: 'N_A',
    });

    const [ordinationData, setOrdinationData] = useState(null);
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
                advancedSettings: {},
                ...parameters
            };

            // Call IPC method to run ordination script
            const result = await window.electron.runOrdinationScript(params)
            
            // Update the ordination data with result
            setOrdinationData(result);
        } catch (error) {
            console.error("Error generating ordination plot:", error);
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
        // Send IPC event to open the ordination PDF
        window.electron.send('open-ordination-location');
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

            <h1>Generate Ordination Plot</h1>

        <div>
            {/* File upload for RData file */}
            <h3>Upload R Data File for Ordination Plot</h3>
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


            {/* Form for Ordination Parameters */}
            <form onSubmit={handleSubmit}>
                {Object.keys(parameters).map((key) => (
                    <div key={key}>
                        <label htmlFor={key}>{key}:</label>
                        {typeof parameters[key] === 'boolean' ? (
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
                    Generate Ordination
                </Button>
            </form> 

            {ordinationData && (
                <div id="ordination=container">
                    <h2>Ordination Plot</h2>
                    <Button
                    variant="contained"
                    color='primary'
                    onClick={handleDownloadClick}
                    >
                        Open Ordination Plot PDF
                    </Button>
                </div>
            )}



        </div>
    );

};

export default Ordination;
