const { contextBridge, ipcRenderer } = require('electron');

contextBridge.exposeInMainWorld('electron', {
  send: (channel, data) => ipcRenderer.send(channel, data),
  on: (channel, func) => ipcRenderer.on(channel, (event, ...args) => func(...args)),
  loadRDataFile: (filePath) => ipcRenderer.invoke('load-rdata-file', filePath),
  runHeatmapScript: (params) => ipcRenderer.invoke('run-heatmap-script', params),
  runOrdinationScript: (params) => ipcRenderer.invoke('run-ordination-script', params),
  runAlphaDiversityScript: (params) => ipcRenderer.invoke('run-alphaDiversity-script', params),
  runRelabundFeaturesScript: (params) => ipcRenderer.invoke('run-relabundFeatures-script', params),
  onParamStr: (callback) => ipcRenderer.on('param-str', (event, paramStr) => callback(paramStr))
});
