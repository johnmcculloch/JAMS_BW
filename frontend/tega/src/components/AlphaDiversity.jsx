import React, { useEffect, useState } from 'react';
import Button from '@mui/material/Button';
import CloudUploadIcon from '@mui/icons-material/CloudUpload';

const AlphaDiversity = () => {
    const [parameters, setParameters] = useState({
        measures: '',
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
    
    const [alphadiversityData, setAlphaDiversityData] = useState(null);
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

            // Call IPC Method to run AlphaDiversity script
            const result = await window.electron.runAlphaDiversityScript(params) // Verify

            // Update the AD data with result
            setAlphaDiversityData(result);
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
            <div style={{ position: 'absolute', top: '10px', right: '10px' }}>
                <button onClick={handleClick}>
                    Go Back to Home Page
                </button>
            </div>

        </div>
    );
};

export default AlphaDiversity;
