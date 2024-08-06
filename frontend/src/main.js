// main.js

const { app, BrowserWindow, ipcMain } = require('electron');
const { exec } = require('child_process');
const path = require('path');

// Handle creating/removing shortcuts on Windows when installing/uninstalling.
if (require('electron-squirrel-startup')) {
  app.quit();
}

//try {
 // require('electron-reloader')(module);
//} catch (_) {}


const createWindow = () => {
  // Create the browser window.
  const mainWindow = new BrowserWindow({
    width: 800,
    height: 600,
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      contextIsolation: true,
      enableRemoteModule: false,
      nodeIntegration: false,
    },
  });

  // and load the index.html of the app.
  mainWindow.loadFile(path.join(__dirname, 'renderer', 'heatmap', 'heatmap.html'));


  // Open the DevTools.
  mainWindow.webContents.openDevTools();
};

// This method will be called when Electron has finished
// initialization and is ready to create browser windows.
app.whenReady().then(() => {
  createWindow();


  app.on('activate', () => {
    if (BrowserWindow.getAllWindows().length === 0) {
      createWindow();
    }
  });
});


app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

// Load r-data-file containing Expvec
ipcMain.handle('load-rdata-file', async (event, filePath) => {
  const script = `
    Rscript -e '
      load("${filePath}");
      list_objs <- ls();
      list_objs <- list_objs[sapply(list_objs, function(x) is.list(get(x)))]
      if (length(list_objs) > 0) {
        first_list <- get(list_objs[1])
        obj_names <- names(first_list)
        writeLines(paste(list_objs[1], obj_names, sep="$"), stdout());
      } else {
        writeLines("", stdout());
      }
    '
  `;

  return new Promise((resolve, reject) => {
    exec(script, (error, stdout, stderr) => {
      if (error) {
        reject(`Error: ${error.message}`);
        return;
      }
      if (stderr) {
        reject(`Stderr: ${stderr}`);
        return;
      }
      const objects = stdout.split('\n').filter(name => name.trim() !== '');
      listObjs = objects
      resolve(objects); // Ensure objects are returned here
    });
  });
});



// get user defined parameters then execute R heatmap script
ipcMain.handle('run-r-script', async (event, params) => {
  // Log the params object for debugging
  console.log('Recieved params:', params);

  const { filePath, ExpObj, ...otherParams } = params;


  // Log individual properties for debugging
  //console.log('filePath:', filePath);
  //console.log('ExpObj:', ExpObj);
  //console.log('otherParams:', otherParams);



  // Construct paramStr dynamically
  const paramStr = `ExpObj = ${ExpObj}, ` +
    Object.entries(otherParams)
      .map(([key, value]) => {
        if (value === "" || value === null || value === 'null' || value === 'NULL') {
          return `${key}=NULL`;
        } else if (typeof value === 'string' && value.startsWith('c(')) { // Check if param is a R variable
          return `${key}=${value}`;
        } else if (typeof value === 'boolean') {
          return `${key}=${value ? 'TRUE' : 'FALSE'}`;
        } else if (typeof value ==='number' || !isNaN(value)) {
          return `${key}=${Number(value)}`
        } else if (typeof value === 'string') {
          return `${key}="${value}"`;
        } else {
          return `${key}=${value}`;
        }
      })

      .join(', ');

  console.log(paramStr);

  // Send paramStr to the renderer process for debugging
  event.sender.send('param-str', paramStr);

  const outputFilePath = path.join(__dirname, 'assets', 'heatmap.pdf');
  const script = `
    Rscript -e '
    suppressPackageStartupMessages({
    load("${filePath}");
    library(JAMS); 
    source("/Users/mossingtonta/Projects/JAMS_BW/R/plot_relabund_heatmap.R"); 
    pdf("${outputFilePath}");
    plot_relabund_heatmap(${paramStr})
    dev.off();
    })'
  `;



  return new Promise((resolve, reject) => {
    exec(script, (error, stdout, stderr) => {
      if (error) {
        reject(`Error: ${error.message}`);
        return;
      }
      if (stderr) {
        reject(`Stderr: ${stderr}`);
        return;
      }
      resolve({ stdout: stdout, imagePath: outputFilePath });
    });
  });
});
