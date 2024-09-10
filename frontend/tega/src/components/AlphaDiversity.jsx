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

    const handleClick = () => {
        window.electron.send('navigate-to', 'home');
    };

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
