const { contextBridge, ipcRenderer } = require('electron');

contextBridge.exposeInMainWorld('electron', {
  send: (channel, data) => ipcRenderer.send(channel, data),
  receive: (channel, func) => ipcRenderer.on(channel, (event, ...args) => func(...args)),
  loadRDataFile: (filePath) => ipcRenderer.invoke('load-rdata-file', filePath),
  invoke: (channel, data) => ipcRenderer.invoke(channel, data),
  runHeatmapScript: (params) => ipcRenderer.invoke('run-heatmap-script', params),
  runOrdinationScript: (params) => ipcRenderer.invoke('run-ordination-script', params),
  runAlphaDiversityScript: (params) => ipcRenderer.invoke('run-alphaDiversity-script', params),
  runRelabundFeaturesScript: (params) => ipcRenderer.invoke('run-relabundFeatures-script', params),
  onParamStr: (callback) => ipcRenderer.on('param-str', (event, paramStr) => callback(paramStr))
});

contextBridge.exposeInMainWorld('updater', {
  onUpdateAvailable: (callback) => {
    console.log('Registering update-available listener');
    ipcRenderer.on('update-available', () => {
      console.log('Received update-available event');
      callback();
    });
  },
  onUpdateDownloaded: (callback) => {
    console.log('Registering update-downloaded listener');
    ipcRenderer.on('update-downloaded', () => {
      console.log('Received update-downloaded event');
      callback();
    });
  },
  installUpdate: () => {
    console.log('Sending install-update command');
    ipcRenderer.send('install-update');
  },
  checkForUpdates: () => ipcRenderer.invoke('check-for-updates'),
  testNotification: (type) => ipcRenderer.invoke('test-update-notification', type)
});

contextBridge.exposeInMainWorld('appInfo', {
  getVersion: () => ipcRenderer.invoke('get-app-version')
});
