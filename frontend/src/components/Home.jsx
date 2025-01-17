import React from 'react';
import { Grid, Paper, Typography, Button } from '@mui/material';

export default function Home({ handleNavigateTo }) {
    return (
        <div>
            {/* Main content area */}
            <Grid container spacing={3}>
                <Grid item xs={12}>
                    <Typography variant="h4" gutterBottom>
                        JAMS Analysis Dashboard
                    </Typography>
                </Grid>

                <Grid item xs={12} sm={6} md={4}>
                    <Paper elevation={3} sx={{ padding: 2 }}>
                        <Typography variant="h6">Heatmap Analysis</Typography>
                        <Button
                            variant="outlined"
                            color="primary"
                            onClick={handleNavigateTo('heatmap')}
                            fullWidth
                        >
                            Go to Heatmap Analysis
                        </Button>
                    </Paper>
                </Grid>

                <Grid item xs={12} sm={6} md={4}>
                    <Paper elevation={3} sx={{ padding: 2 }}>
                        <Typography variant="h6">Ordination Analysis</Typography>
                        <Button
                            variant="outlined"
                            color="primary"
                            onClick={handleNavigateTo('ordination')}
                            fullWidth
                        >
                            Go to Ordination Analysis
                        </Button>
                    </Paper>
                </Grid>

                <Grid item xs={12} sm={6} md={4}>
                    <Paper elevation={3} sx={{ padding: 2 }}>
                        <Typography variant="h6">Alpha Diversity Analysis</Typography>
                        <Button
                            variant="outlined"
                            color="primary"
                            onClick={handleNavigateTo('alphadiversity')}
                            fullWidth
                        >
                            Go to Alpha Diversity Analysis
                        </Button>
                    </Paper>
                </Grid>

                <Grid item xs={12} sm={6} md={4}>
                    <Paper elevation={3} sx={{ padding: 2 }}>
                        <Typography variant="h6">Relabund Features Analysis</Typography>
                        <Button
                            variant="outlined"
                            color="primary"
                            onClick={handleNavigateTo('relabundfeatures')}
                            fullWidth
                        >
                            Go to Relabund Features Analysis
                        </Button>
                    </Paper>
                </Grid>
            </Grid>
        </div>
    );
}
