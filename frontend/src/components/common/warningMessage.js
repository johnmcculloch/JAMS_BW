import React, { useState } from 'react';
import { Alert, Snackbar } from '@mui/material';

const WarningSnackbar = ({ open, message, onClose }) => {
    return (
        <Snackbar
            open={open}
            autoHideDuration={6000}
            onClose={onClose}
            anchorOrigin={{ vertical: 'top', horizontal: 'center' }}
        >
            <Alert
                onClose={onClose}
                severity='warning'
                sx={{ width: '100%' }}
            >
                {message}
            </Alert>
        </Snackbar>  
    );
};

// Validation functions
export const validateRequiredFields = (fields) => {
    // Check each field in the object
    for (const [fieldName, value] of Object.entries(fields)) {
        if (!value || (Array.isArray(value) && value.length === 0)) {
            return {
                isValid: false,
                message: `Please provide a ${fieldName.replace(/([A-Z])/g, ' $1').toLowerCase()}`,
                field: fieldName
            };
        }
    }

    return {
        isValid: true
    };
};

// Custom hook for warning management
export const useWarningState = () => {
    const [showWarning, setShowWarning] = useState(false);
    const [warningMessage, setWarningMessage] = useState('');

    const handleCloseWarning = () => {
        setShowWarning(false);
    };

    const showWarningMessage = (message) => {
        setWarningMessage(message);
        setShowWarning(true);
    };

    return {
        showWarning,
        warningMessage,
        handleCloseWarning,
        showWarningMessage
    };
};

export default WarningSnackbar;
