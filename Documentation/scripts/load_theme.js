function load_theme() {
  try {
    button = document.getElementById("Theme_button_logo");
    if (localStorage.getItem("theme") == "light") {
      document.getElementById('css_style').href = 'styles/light.css';
      var images = document.getElementsByTagName('img'); 
      for(var i = 0; i < images.length; i++) {
        images[i].style.filter = "invert(0)";
      }
      button.src = "scripts/Font-awesome/svgs/solid/moon.svg";
    }
    else {
      document.getElementById('css_style').href = 'styles/dark.css';
      var images = document.getElementsByTagName('img'); 
      for(var i = 0; i < images.length; i++) {
        if (images[i].id == "logoIIT")
        images[i].style.filter = "invert(1) hue-rotate(180deg) brightness(1.5)";
        else {
          if (images[i].id == "Logo_VFlexP")
            images[i].style.filter = "invert(0)";
          else
            images[i].style.filter = "invert(1) hue-rotate(180deg) brightness(1.5)";
        }
      } 
      button.src = "scripts/Font-awesome/svgs/solid/sun.svg";
    }
  }
  catch(err) {
    console.log(err);
  }
}

function load_theme_pre() {
  if (localStorage.getItem("theme") == "dark") {
    document.body.style.backgroundColor = "black";
  }
  else {
    document.body.style.backgroundColor = "white";
  }
}