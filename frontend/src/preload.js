const { contextBridge, ipcRenderer } = require('electron');

contextBridge.exposeInMainWorld('electron', {
  send: (channel, data) => ipcRenderer.send(channel, data),
  on: (channel, func) => ipcRenderer.on(channel, (event, ...args) => func(...args)),
  loadRDataFile: (filePath) => ipcRenderer.invoke('load-rdata-file', filePath),
  runRScript: (params) => ipcRenderer.invoke('run-r-script', params),
  onParamStr: (callback) => ipcRenderer.on('param-str', (event, paramStr) => callback(paramStr))
});
