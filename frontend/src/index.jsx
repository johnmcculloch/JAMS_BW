import * as React from 'react';
import { useState, useEffect } from 'react';
import { createRoot } from 'react-dom/client';
import { ThemeProvider, createTheme, CssBaseline, useMediaQuery, Drawer, List, ListItem, ListItemText, IconButton, Box, Menu, MenuItem, Switch } from '@mui/material';
import MenuIcon from '@mui/icons-material/Menu';
import MoreVertIcon from '@mui/icons-material/MoreVert';
import ExpandLess from '@mui/icons-material/ExpandLess';
import ExpandMore from '@mui/icons-material/ExpandMore';
import Home from './components/Home.jsx';
import Heatmap from './components/Heatmap.jsx';
//import Ordination from './components/Ordination.jsx';
//import AlphaDiversity from './components/AlphaDiversity.jsx';
//import RelabundFeatures from './components/RelabundFeatures.jsx';

function App() {
    const [currentPage, setCurrentPage] = useState('home');
    const [isDrawerOpen, setIsDrawerOpen] = useState(false);
    const [isSettingsOpen, setIsSettingsOpen] = useState(false);
    const [anchorEl, setAnchorEl] = useState(null);
    const [darkMode, setDarkMode] = useState(false);

    const prefersDarkMode = useMediaQuery('(prefers-color-scheme: dark)');

    useEffect(() => {
        setDarkMode(prefersDarkMode);
    }, [prefersDarkMode]);

    const theme = createTheme({
        palette: {
            mode: darkMode ? 'dark' : 'light',
        },
    });

    const toggleDrawer = (open) => (event) => {
        if (event.type === 'keydown' && (event.key === 'Tab' || event.key === 'Shift')) {
            return;
        }
        setIsDrawerOpen(open);
    };

    const handleNavigateTo = (page) => () => {
        setCurrentPage(page);
        window.electron.send('navigate-to', page);
    };

    const handleToggleSettings = () => {
        setIsSettingsOpen(!isSettingsOpen);
    };

    const handleMenuOpen = (event) => {
        setAnchorEl(event.currentTarget);
    };

    const handleMenuClose = () => {
        setAnchorEl(null);
    };

    const handleToggleDarkMode = () => {
        setDarkMode((prevMode) => !prevMode);
    };

    const renderNavigationList = () => (
        <Box sx={{ width: 250 }} role="presentation" onClick={toggleDrawer(false)} onKeyDown={toggleDrawer(false)}>
            <List>
                <ListItem button onClick={handleNavigateTo('home')}>
                    <ListItemText primary="Home" />
                </ListItem>
                <ListItem button onClick={handleNavigateTo('heatmap')}>
                    <ListItemText primary="Heatmap Analysis" />
                </ListItem>
                <ListItem button onClick={handleNavigateTo('ordination')}>
                    <ListItemText primary="Ordination Analysis" />
                </ListItem>
                <ListItem button onClick={handleNavigateTo('alphadiversity')}>
                    <ListItemText primary="Alpha Diversity Analysis" />
                </ListItem>
                <ListItem button onClick={handleNavigateTo('relabundfeatures')}>
                    <ListItemText primary="Relabund Features Analysis" />
                </ListItem>
                <ListItem button onClick={handleToggleSettings}>
                    <ListItemText primary="Settings" />
                    {isSettingsOpen ? <ExpandLess /> : <ExpandMore />}
                </ListItem>
                {isSettingsOpen && (
                    <List component="div" disablePadding>
                        <ListItem button sx={{ pl: 4 }}>
                            <ListItemText primary="Profile Settings" />
                        </ListItem>
                        <ListItem button sx={{ pl: 4 }}>
                            <ListItemText primary="App Settings" />
                        </ListItem>
                    </List>
                )}
            </List>
        </Box>
    );

    return (
        <ThemeProvider theme={theme}>
            <CssBaseline />
            <div>
                {/* Sidebar menu button */}
                <IconButton onClick={toggleDrawer(true)} sx={{ position: 'absolute', top: 16, left: 16 }}>
                    <MenuIcon />
                </IconButton>

                {/* Drawer for navigation sidebar */}
                <Drawer open={isDrawerOpen} onClose={toggleDrawer(false)}>
                    {renderNavigationList()}
                </Drawer>

                {/* Main content area */}
                <div style={{ marginLeft: 250, padding: 16 }}>
                    {currentPage === 'home' && <Home handleNavigateTo={handleNavigateTo} />}
                    {currentPage === 'heatmap' && <Heatmap handleNavigateTo={handleNavigateTo} />}
                    {currentPage === 'ordination' && <Ordination handleNavigateTo={handleNavigateTo} />}
                    {currentPage === 'alphadiversity' && <AlphaDiversity handleNavigateTo={handleNavigateTo} />}
                    {currentPage === 'relabundfeatures' && <RelabundFeatures handleNavigateTo={handleNavigateTo} />}
                </div>

                {/* Settings Icon in top-right corner */}
                <Box sx={{ position: 'absolute', top: 16, right: 16 }}>
                    <IconButton
                        aria-label='settings'
                        aria-controls='settings-menu'
                        aria-haspopup='true'
                        onClick={handleMenuOpen}
                    >
                        <MoreVertIcon />
                    </IconButton>

                    {/* Settings Menu */}
                    <Menu
                        id='settings-menu'
                        anchorEl={anchorEl}
                        open={Boolean(anchorEl)}
                        onClose={handleMenuClose}
                    >
                        <MenuItem>
                            <ListItemText primary="Dark Mode" />
                            <Switch
                                checked={darkMode}
                                onChange={handleToggleDarkMode}
                                inputProps={{ 'aria-label': 'controlled' }}
                            />
                        </MenuItem>
                    </Menu>
                </Box>
            </div>
        </ThemeProvider>
    );
}

const root = createRoot(document.getElementById('root'));
root.render(<App />);
