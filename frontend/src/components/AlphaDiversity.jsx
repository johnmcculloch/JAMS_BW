import React, { useEffect, useState, useRef } from 'react';
import Button from '@mui/material/Button';
import CloudUploadIcon from '@mui/icons-material/CloudUpload';
import Box from '@mui/material/Box';
import useAutoScroll from './common/useAutoScroll';
import LoadingIndicator from './common/LoadingIndicator';
import WarningSnackbar, { validateRequiredFields, useWarningState } from './common/warningMessage';

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
    const [loading, setLoading] = useState(false);
    const [generatingAD, setGeneratingAD] = useState(false);

    const resultRef = useRef(null);
    const loadingRef = useRef(null);

    // auto scroll to results when ready
    useAutoScroll(alphaDiversityData, resultRef);

    // auto scroll to loading indicator
    useAutoScroll(generatingAD, loadingRef);

    // use warning state hook
    const { showWarning, warningMessage, handleCloseWarning, showWarningMessage } = useWarningState();

    const handleChange = (e) => {
        const { name, value, type, checked } = e.target;
        setParameters({
            ...parameters,
            [name]: type === 'checkbox' ? checked : value
        });
    };

    const handleSubmit = async (e) => {
        e.preventDefault();

        const validationResult = validateRequiredFields({
            'R data file': filePath,
            'Summarized Experiment Object': selectedObj
        });

        if (!validationResult.isValid) {
            showWarningMessage(validationResult.message);
            return;
        }

        console.log('Selected ExpObj:', selectedObj);
        try {
            setLoading(true);
            setGeneratingAD(true);
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
            setWarningMessage(`Error: ${error.message || 'Failed to generate alpha diversity plot'}`);
            setShowWarning(true);
        } finally {
            setLoading(false);
            setGeneratingAD(false);
        }
    };

    const handleFileUpload = async () => {
        try {
            setLoading(true);
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
        } finally {
            setLoading(false);
        }
    };

    const handleObjSelect = (e) => {
        const value = e.target.value;
        setSelectedObj(value);
        console.log('Selected object state:', value);
    };

    const handleDownloadClick = () => {
        // Send IPC event to open the alphadiversity PDF
        window.electron.send('open-alphadiversity-location');
    };

    useEffect (() => {
        // Set the default selection if objects array changes
        if (objects.length > 0 && !selectedObj) {
            setSelectedObj(objects[0]);
        }
    }, [objects, selectedObj]);

    return (
        <div>
            <Button
                onClick={handleNavigateTo('home')}
                sx={{
                    position: 'absolute',
                    top: '60px',
                    right: '10px',
                }}
            >
                Go Back to Home Page
            </Button>

            {/* warning message */}
            <WarningSnackbar
                open={showWarning}
                message={warningMessage}
                onClose={handleCloseWarning}
            />

            <h1>Generate AlphaDiversity Plot</h1>
        <div>
            {/* File upload for RData file */}
            <h3>Upload RData File for AlphaDiversity Plot</h3>
            <Button
                component="label"
                variant="contained"
                startIcon={<CloudUploadIcon />}
                onClick={handleFileUpload}
            >
                Upload RData File
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

        {/* Form for AlphaDiversity Parameters */}
        <form onSubmit={handleSubmit}>
            {Object.keys(parameters).map((key) => (
                <div key={key}>
                    <label htmlFor={key}>{displayNames[key] || key}:</label>
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

                    ): typeof parameters[key] === 'boolean' ? (
                        <input
                        type='checkbox'
                        id={key}
                        name={key}
                        checked={parameters[key]}
                        onChange={handleChange}
                        />
                    ): (
                        <input
                        type='text'
                        id={key}
                        name={key}
                        value={parameters[key]}
                        onChange={handleChange}
                        />
                    )}
                </div>
            ))}
            <Button type='submit' variant='outlined' color='primary'>
                Generate Alpha Diversity Plot
            </Button>
        </form>

        {/* Add loading indicator */}
        {loading && <LoadingIndicator ref={generatingAD ? loadingRef : null} />}

        {alphaDiversityData && !loading && (
            <div id='alphadiversity=container' ref={resultRef}>
                <h2>AlphaDiversity Plot</h2>
                <Button
                variant='contained'
                color='primary'
                onClick={handleDownloadClick}
                >
                    Open AlphaDiversity Plot PDF
                </Button>
            </div>
        )}
    </div>
    );
};

export default AlphaDiversity;
