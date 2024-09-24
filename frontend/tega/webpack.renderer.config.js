const path = require('path');
const CopyWebpackPlugin = require('copy-webpack-plugin');
const { DefinePlugin } = require('webpack'); // Import DefinePlugin
const rules = require('./webpack.rules');

rules.push({
  test: /\.css$/,
  use: [{ loader: 'style-loader' }, { loader: 'css-loader' }],
});

module.exports = {
  entry: './src/index.jsx',
  output: {
    filename: 'renderer.bundle.js',
    path: path.resolve(__dirname, 'dist'),
  },
  module: {
    rules,
  },
  plugins: [
    new CopyWebpackPlugin({
      patterns: [
        { from: path.resolve(__dirname, 'src/index.html'), to: path.resolve(__dirname, 'dist') },
      ],
    }),
    new DefinePlugin({
      __dirname: JSON.stringify(path.resolve(__dirname, 'dist')), // Define __dirname
    }),
  ],
  target: 'electron-renderer',
  resolve: {
    extensions: ['.js', '.jsx'],
  },
};
