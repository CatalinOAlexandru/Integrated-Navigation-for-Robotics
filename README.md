# Integrated-Navigation-for-Robotics

In this report, we document our solution for the coursework. We are looking at the problem of integrated navigation of a robotic lawnmower, using GNSS and Dead Reckoning (DR) sensor data. Our approach has 4 pillars: 1) initialising the position using GNSS 2) Outlier Detection in the data 3) use Kalman filter for the GNSS only positioning 4) use DR data to create a DR only solution, then integrate it with the GNSS positions using Kalman filter.

# Some Results:

GNSS Only
![GNSS_only](https://user-images.githubusercontent.com/33178694/121712695-ed1dde80-cad3-11eb-99fb-cd03c425aee8.png)

Dead reckoning
![DR_only](https://user-images.githubusercontent.com/33178694/121712723-f6a74680-cad3-11eb-86ea-4dd3304b8549.png)

Integrated with Outlier Detection
![Integraded_only_final](https://user-images.githubusercontent.com/33178694/121712851-176f9c00-cad4-11eb-8da7-e4ba5d83ef2e.png)
