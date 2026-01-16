Project in progress.

Reduce the wave-like effect near edges in the image, mainly due to the distortion of light crossing through the air with various temperatures. It often happens in motorsport photography and plane spotting. 

The program automatically detects the edge in selected area (make sure only include one piece of edge in every brush). It then calculates the best-fit line (with DOF=3) and render the original edge into the refined edge.

Run the project (I run it in Visual Studio). Use the mouse to do the following operations:

Zoom in/out: scroll the mouse (do not scroll too fast. Seeing bugs when keep zooming out when reaching the limit)

Move: right button down, then move the mouse.

Assign edge to be edited: Left button down, then move along the edge of the image.

Undo: mid button (wheel) down.

![currentimg](https://github.com/user-attachments/assets/015cfbd0-00d5-4ff2-b48f-e691c330dce3)
![rrsult](https://github.com/user-attachments/assets/b116dbe1-14a9-4cb0-8bd6-5eb28b89b848)
