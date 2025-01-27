const { app, BrowserWindow, ipcMain, dialog, shell } = require('electron');
const path = require('path');
const { exec } = require('child_process');
const fs = require('fs');

let mainWindow;

console.log(`app.isPackaged: ${app.isPackaged}`); // debugging
let isDev = !app.isPackaged;

console.log(`isDev: ${isDev}`);

function createWindow() {
  mainWindow = new BrowserWindow({
    width: 800,
    height: 600,
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
    },
  });

  const startUrl = isDev
    ? 'http://localhost:1234'
    : `file://${path.join(__dirname, 'dist', 'index.html')}`;
  console.log(`Loading URL: ${startUrl}`); 
  mainWindow.loadURL(startUrl);

  mainWindow.on('closed', function () {
    mainWindow = null;
  });
}

app.on('ready', createWindow);

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

app.on('activate', () => {
  if (mainWindow === null) {
    createWindow();
  }
});

// ** Main Processes ** 

// Handle Page Navigation IPC event
ipcMain.on('navigate-to', (event, page) => {
  // Send the navigation event to the renderer
  mainWindow.webContents.send('navigate-to', page);
});


// Function to handle RData session files
ipcMain.handle('open-file-dialog', async () => {
  const { canceled, filePaths } = await dialog.showOpenDialog({
    properties: ['openFile'],
    filters: [
      { name: 'RData Files', extensions: ['rdata', 'rda'] }
    ]
  });

  if (canceled) {
    return null;
  } else {
    return filePaths[0];
  }
});


// Load user provided r-data-file containing Expvec
ipcMain.handle('load-rdata-file', async (event, filePath) => {
  return new Promise((resolve, reject) => {
    const rscriptPath = '/usr/local/bin/Rscript'; // locate R on user's system
    const command = `
    ${rscriptPath} -e '
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
  exec(command, (error, stdout, stderr) => {
    if (error) {
      console.error(`exec error: ${error}`);
      reject(`Error: ${error.message}`);
      return;
    }
    if (stderr) {
      console.error(`stderr: ${stderr}`);
      return;
    }
    const objects = stdout.split('\n').filter(name => name.trim() !== '');
    //listObjs = objects (VERIFY REDUNDANCY)
    resolve(objects); // Ensure objects are returned here
    });
  });
});


// HEATMAP Script 
ipcMain.handle('run-heatmap-script', async (event, params) => {
  // Log the params object for debugging
  console.log('Recieved params:', params);

  const { filePath, ExpObj, advancedSettings, ...otherParams } = params;

  // Construct paramStr dynamically to account for anything the user inputs
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

// Run plot_relabund_heatmap command with user-defined parameters
  const outputDir = path.join(app.getPath('userData'), 'assets');
  const outputFilePath = path.join(outputDir, 'heatmap.pdf');
  const rscriptPath = '/usr/local/bin/Rscript';

  // Create output directory and log results
  try {
    if (!fs.existsSync(outputDir)) {
      fs.mkdirSync(outputDir, { recursive: true });
      console.log(`Created directory: ${outputDir}`);
    }
    console.log(`Output directory exists: ${outputDir}`);
  } catch (err) {
    console.error(`Error creating directory: ${err}`);
  }

  // Escape spaces in file path for R
  const escapedOutputPath = outputFilePath.replace(/ /g, '\\ ');

  const script = `
   ${rscriptPath} -e '
    suppressPackageStartupMessages({
    load("${filePath}");
    library(JAMS); 
    source("/Users/mossingtonta/Projects/JAMS_BW_DEV/R/plot_relabund_heatmap.R"); 
    pdf("${escapedOutputPath}");
    plot_relabund_heatmap(${paramStr})
    dev.off();
    })'
  `;

  console.log('Executing command:', script);

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

// IPC handler for opening the heatmap file
ipcMain.on('open-heatmap-location', (event, filePath) => {
  const outputFilePath = path.join(app.getPath('userData'), 'assets', 'heatmap.pdf');
  shell.openPath(outputFilePath);
});


// get user defined parameters then execute R ordination plot script
ipcMain.handle('run-ordination-script', async (event, params) => {
  // Log the params object for debugging
  console.log('Recieved params:', params);

  const { filePath, ExpObj, ...otherParams } = params;

  // Construct paramStr dynamically to account for anything the user inputs
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

// Run plot_Ordination command with user-defined parameters
  const outputFilePath = path.join(__dirname, 'assets', 'ordination.pdf');
  const script = `
    Rscript -e '
    suppressPackageStartupMessages({
    load("${filePath}");
    library(JAMS); 
    source("/Users/mossingtonta/Projects/JAMS_BW_DEV/R/plot_Ordination.R"); 
    pdf("${outputFilePath}", paper = "a4r");
    print(plot_Ordination(${paramStr}))
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


// IPC handler for opening the ordination plot file (ordination.js)
ipcMain.on('open-ordination-location', (event) => {
  const outputFilePath = path.join(__dirname, 'assets', 'ordination.pdf');
  shell.openPath(outputFilePath).then((error) => {
    if (error) {
      console.error('Failed to open file location:', error);
    } else {
      console.log('File location opened successfully!');
    }
  });
});


// get user defined parameters then execute R alpha diversity plot script
ipcMain.handle('run-alphaDiversity-script', async (event, params) => {
  // Log the params object for debugging
  console.log('Recieved params:', params);

  const { filePath, ExpObj, ...otherParams } = params;

  // Construct paramStr dynamically to account for anything the user inputs
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

// Run plot_alphadiversity command with user-defined parameters
  const outputFilePath = path.join(__dirname, 'assets', 'alphaDiversity.pdf');
  const script = `
    Rscript -e '
    suppressPackageStartupMessages({
    suppressWarnings({
      load("${filePath}");
      library(JAMS); 
      source("/Users/mossingtonta/Projects/JAMS_BW_DEV/R/plot_alpha_diversity.R"); 
      pdf("${outputFilePath}", paper = "a4r");
      print(plot_alpha_diversity(${paramStr}))
      dev.off();
      })
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

// IPC handler for opening the alpha diversity file (alpha_diversity.js)
ipcMain.on('open-alphadiversity-location', (event) => {
  const outputFilePath = path.join(__dirname, 'assets', 'alphaDiversity.pdf');
  shell.openPath(outputFilePath).then((error) => {
    if (error) {
      console.error('Failed to open file location:', error);
    } else {
      console.log('File location opened successfully!');
    }
  });
});

// get user defined parameters then execute R relabund features plot script
ipcMain.handle('run-relabundFeatures-script', async (event, params) => {
  // Log the params object for debugging
  console.log('Recieved params:', params);

  const { filePath, ExpObj, ...otherParams } = params;

  // Construct paramStr dynamically to account for anything the user inputs
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

// Run plot_Ordination command with user-defined parameters
  const outputFilePath = path.join(__dirname, 'assets', 'relabundFeatures.pdf');
  const script = `
    Rscript -e '
    suppressPackageStartupMessages({
    suppressWarnings({
      load("${filePath}");
      library(JAMS);  
      pdf("${outputFilePath}", paper = "a4r");
      print(plot_relabund_features(${paramStr}))
      dev.off();
      })
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


// IPC handler for opening the relabund features file (plot_relabund_features.js)
ipcMain.on('open-RelabundFeatures-location', (event) => {
  const outputFilePath = path.join(__dirname, 'assets', 'relabundFeatures.pdf');
  shell.openPath(outputFilePath).then((error) => {
    if (error) {
      console.error('Failed to open file location:', error);
    } else {
      console.log('File location opened successfully!');
    }
  });
});
