# EBWFNet
Convolutional neural network for improved event-based Shack-Hartmann wavefront reconstruction

Processing pipeline:
1) Convert raw event data to .mat files (save_events_as_mat.m)
2) Identify hot pixels (identify_hot_pixels.m)
3) Select valid subapertures (match_spots.m)
4) Process frame data (process_frame.m)
5) Process event data (process_event.m)
6) Format data (CNN_FormatData.m)
7) Train CNN (train_CNN_CustomLoop.m)
