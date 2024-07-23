const { contextBridge, ipcRenderer } = require('electron');

contextBridge.exposeInMainWorld('electron', {
  loadRDataFile: (filePath) => ipcRenderer.invoke('load-rdata-file', filePath),
  runRScript: (params) => ipcRenderer.invoke('run-r-script', params)
});
