import cv2

# Global variable to store clicked coordinates
coordinates = []

# Mouse callback function to record coordinates
def click_event(event, x, y, flags, param):
    if event == cv2.EVENT_LBUTTONDOWN:
        coordinates.append((x, y))
        print(f"Clicked at: ({x}, {y})")

# Function to process each snapshot
def process_snapshot(image_path):
    # Load the image
    image = cv2.imread(image_path)

    # Resize the image to half its original size
    height, width = image.shape[:2]
    resized_image = cv2.resize(image, (width // 2, height // 2))

    # Display the resized image
    cv2.imshow(f"Resized {image_path}", resized_image)
    cv2.setMouseCallback(f"Resized {image_path}", click_event)

    # Wait until a key is pressed
    cv2.waitKey(0)
    cv2.destroyAllWindows()

    # Optionally, save the coordinates to a file
    with open(f'coordinates_{image_path}.txt', 'w') as file:
        for coord in coordinates:
            file.write(f"{coord[0]}, {coord[1]}\n")

# Process all snapshots
for i in range(4):  # snapshot0.jpg to snapshot3.jpg
    image_path = f"snapshot{i:01d}.jpg"
    coordinates.clear()  # Clear the coordinates list for the next image
    process_snapshot(image_path)
