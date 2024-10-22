const path = require('path');

module.exports = {
  // This is the main entry point for your application, it's the first file that runs in the main process.
  entry: {
    main: './src/main.js',
    preload: './src/preload.js',
  },
  output: {
    filename: '[name].bundle.js',
    path: path.resolve(__dirname, 'dist'),
  },
  module: {
    rules: require('./webpack.rules'),
  },
};
