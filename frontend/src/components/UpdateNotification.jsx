import React, { useEffect, useState } from 'react';
import { Alert, Button, Snackbar } from '@mui/material';

const UpdateNotification = () => {
  const [updateAvailable, setUpdateAvailable] = useState(false);
  const [updateDownloaded, setUpdateDownloaded] = useState(false);

  useEffect(() => {
    // Listen for update events from the main process
    if (window.updater) {
      window.updater.onUpdateAvailable(() => {
        console.log('Update available');
        setUpdateAvailable(true);
      });

      window.updater.onUpdateDownloaded(() => {
        console.log('Update downloaded');
        setUpdateAvailable(false);
        setUpdateDownloaded(true);
      });
    }
  }, []);

  const handleInstallUpdate = () => {
    window.updater.installUpdate();
  };

  return (
    <>
      <Snackbar 
        open={updateAvailable} 
        anchorOrigin={{ vertical: 'bottom', horizontal: 'center' }}
      >
        <Alert severity="info">
          A new version is being downloaded...
        </Alert>
      </Snackbar>
      
      <Snackbar 
        open={updateDownloaded} 
        anchorOrigin={{ vertical: 'bottom', horizontal: 'center' }}
      >
        <Alert 
          severity="success" 
          action={
            <Button color="inherit" size="small" onClick={handleInstallUpdate}>
              INSTALL NOW
            </Button>
          }
        >
          Update ready to install
        </Alert>
      </Snackbar>
    </>
  );
};

export default UpdateNotification;
