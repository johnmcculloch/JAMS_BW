// src/components/VersionDisplay.jsx
import React, { useState, useEffect } from 'react';
import { Box, Typography } from '@mui/material';

const VersionDisplay = () => {
  const [version, setVersion] = useState('');

  useEffect(() => {
    if (window.appInfo) {
      setVersion(window.appInfo.getVersion());
    }
  }, []);

  return (
    <Box
      sx={{
        position: 'absolute',
        bottom: 8,
        right: 16,
        opacity: 0.7,
      }}
    >
      <Typography variant="caption" color="textSecondary">
        Version {version}
      </Typography>
    </Box>
  );
};

export default VersionDisplay;
