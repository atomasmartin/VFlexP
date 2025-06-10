// Wait for the page to load
document.addEventListener('DOMContentLoaded', function() {

  const mainContent = document.getElementById('main_content');
  // Get the sidebar element
  var sidebar = document.getElementById('sidebar');

  // Create the hamburger button
  const hamburgerButton = document.createElement('button');
  hamburgerButton.classList.add('hamburger-icon');
  hamburgerButton.innerHTML = '&#9776;'; // Unicode character for hamburger icon

  // Insert the button as the first child of main_content
  if (mainContent) {
    mainContent.prepend(hamburgerButton);
  }

  // Add event listener to toggle the sidebar
  if (hamburgerButton && sidebar) {
    hamburgerButton.addEventListener('click', function() {
      sidebar.classList.toggle('sidebar-visible');
    });

    // Optional: Close sidebar if clicking outside of it on small screens
    document.addEventListener('click', function(event) {
      const isClickInsideSidebar = sidebar.contains(event.target);
      const isClickOnHamburger = hamburgerButton.contains(event.target);
      const isSidebarVisible = sidebar.classList.contains('sidebar-visible');

      if (isSidebarVisible && !isClickInsideSidebar && !isClickOnHamburger) {
        sidebar.classList.remove('sidebar-visible');
      }
    });
  }

    // Get the right sidebar element
  var right_sidebar = document.getElementById('right_sidebar');

  // Create the hamburger button
  const right_hamburgerButton = document.createElement('button');
  right_hamburgerButton.id = 'right_hamburgerButton';
  right_hamburgerButton.classList.add('hamburger-icon');
  right_hamburgerButton.innerHTML = '&#9776;'; // Unicode character for hamburger icon

  // Insert the button as the first child of main_content
  if (mainContent) {
    mainContent.prepend(right_hamburgerButton);
  }

  // Add event listener to toggle the sidebar
  if (right_hamburgerButton && right_sidebar) {
    right_hamburgerButton.addEventListener('click', function() {
      right_sidebar.classList.toggle('right_sidebar-visible');
    });

    // Optional: Close right_sidebar if clicking outside of it on small screens
    document.addEventListener('click', function(event) {
      const isClickInsideSidebar = right_sidebar.contains(event.target);
      const isClickOnHamburger = right_hamburgerButton.contains(event.target);
      const isSidebarVisible = right_sidebar.classList.contains('right_sidebar-visible');

      if (isSidebarVisible && !isClickInsideSidebar && !isClickOnHamburger) {
        right_sidebar.classList.remove('right_sidebar-visible');
      }
    });
  }


  // Get all the headers in the document
  var headers = document.querySelectorAll('h1, h2, h3');
  
  // Space
  var space = document.createElement('br');
  sidebar.appendChild(space);
  
  // Create Home button
  var home = document.createElement('button');
  home.id = 'Home_button';
  home.addEventListener('click', function() {
    location.href='doc_Home.html';
  });
  var home_icon = document.createElement('img');
  home_icon.src = "scripts/Font-awesome/svgs/solid/house.svg";
  home_icon.style.width = '1cm';
  home.appendChild(home_icon);
  sidebar.appendChild(home);
  
  // Horizontal space between buttons
  var space = document.createTextNode('\u00A0\u00A0');
  sidebar.appendChild(space);
  
  // Create theme button
  var theme = document.createElement('button');
  theme.id = 'Theme_button';
  theme.addEventListener('click', function() {
    switch_themes();
    try {
      requestAnimationFrame(draw);
    } catch {
      console.log('No canvas to draw');
    }
  });
  var theme_icon = document.createElement('img');
  theme_icon.id = 'Theme_button_logo';
  theme_icon.style.width = '1cm';
  theme_icon.style.height = '0.7cm';
  theme.appendChild(theme_icon);
  sidebar.appendChild(theme);
  
  // Space
  var space = document.createElement('br');
  sidebar.appendChild(space);
  
  // Space
  var space = document.createElement('br');
  sidebar.appendChild(space);
  
  // Create an image element for the icon
  var icon = document.createElement('img');
  icon.src = 'figures/IIT.png';
  icon.alt = 'IIT logo';
  icon.style.width = '7cm';
  icon.id = 'logoIIT';
  
  // Create an anchor element for the link
  var link = document.createElement('a');
  link.href = 'https://www.iit.comillas.edu';
  
  // Append the image to the link
  link.appendChild(icon);
  
  // Append the link to the sidebar
  sidebar.appendChild(link);
  
  // Navigation title
  var nav = document.createElement('p');
  nav.className = 'nolink_sidebar';
  nav.textContent = 'Page navigation';
  sidebar.appendChild(nav);
  
  // Loop through each header
  headers.forEach(function(header) {
    // Create a new sidebar element for the header
    var sidebarItem = document.createElement('p');
    
    // Create a span element for the header text
    var headerText = document.createElement('span');
    headerText.textContent = header.textContent;
    
    // Append the header text to the sidebar element
    sidebarItem.appendChild(headerText);
    
    // Set the class name of the sidebar element based on the header level
    switch (header.tagName.toLowerCase()) {
      case 'h1':
      sidebarItem.className = 'sidebar-h1';
      break;
      case 'h2':
      sidebarItem.className = 'sidebar-h2';
      break;
      case 'h3':
      sidebarItem.className = 'sidebar-h3';
      break;
      default:
      sidebarItem.className = 'sidebar-item';
    }
    
    // Add a click event listener to scroll to the corresponding header
    sidebarItem.addEventListener('click', function() {
      header.scrollIntoView({ behavior: 'smooth' });
    });
    
    // Set the indentation based on the header level
    if (header.tagName.toLowerCase() === 'h2') {
      sidebarItem.style.marginLeft = '1cm'; // Adjust the indentation as needed
    } else if (header.tagName.toLowerCase() === 'h3') {
      sidebarItem.style.marginLeft = '2cm'; // Adjust the indentation as needed
    }
    
    // Append the sidebar element to the sidebar
    sidebar.appendChild(sidebarItem);
  });
  
  // Now for the right sidebar
  var right_sidebar = document.getElementById('right_sidebar');
  
  // Space
  var space = document.createElement('br');
  right_sidebar.appendChild(space);
  var space = document.createElement('br');
  right_sidebar.appendChild(space);
  var space = document.createElement('br');
  right_sidebar.appendChild(space);
  
  // Create an image element for the icon
  var icon = document.createElement('img');
  icon.src = 'figures/Logo.svg';
  icon.style.width = '6cm';
  icon.id = 'Logo_VFlexP';
  
  right_sidebar.appendChild(icon);
  
  // Functionalities of the tool
  var Functionalities = document.createElement('p');
  Functionalities.textContent = 'Functionalities of the tool';
  Functionalities.className = 'nolink_sidebar';
  right_sidebar.appendChild(Functionalities);

  var Tools = document.createElement('p');
  Tools.textContent = 'Basics';
  Tools.className = 'nolink_sidebar';
  right_sidebar.appendChild(Tools);
  
  var pdqtransform = document.createElement('p');
  pdqtransform.style.marginLeft = '1cm';
  var pdqtransform_link = document.createElement('a');
  pdqtransform_link.href = 'doc_dqtransformation.html';
  pdqtransform_link.textContent = 'd-q transformation';
  pdqtransform.appendChild(pdqtransform_link);
  right_sidebar.appendChild(pdqtransform);
  
  var gstruct = document.createElement('p');
  gstruct.style.marginLeft = '1cm';
  var gstruct_link = document.createElement('a');
  gstruct_link.href = 'doc_grid_structure.html';
  gstruct_link.textContent = 'Grid structure';
  gstruct.appendChild(gstruct_link);
  right_sidebar.appendChild(gstruct);

  var load_and_save = document.createElement('p');
  load_and_save.style.marginLeft = '1cm';
  var load_and_save_link = document.createElement('a');
  load_and_save_link.href = 'doc_load_and_save.html';
  load_and_save_link.textContent = 'Load and save';
  load_and_save.appendChild(load_and_save_link);
  right_sidebar.appendChild(load_and_save);

  var Static = document.createElement('p');
  Static.textContent = 'Static';
  Static.className = 'nolink_sidebar';
  right_sidebar.appendChild(Static);

  var powerflow = document.createElement('p');
  powerflow.style.marginLeft = '1cm';
  var powerflow_link = document.createElement('a');
  powerflow_link.href = 'doc_powerflow.html';
  powerflow_link.textContent = 'Power flow';
  powerflow.appendChild(powerflow_link);
  right_sidebar.appendChild(powerflow);

  var initialise_units = document.createElement('p');
  initialise_units.style.marginLeft = '1cm';
  var initialise_units_link = document.createElement('a');
  initialise_units_link.href = 'doc_initialise_units.html';
  initialise_units_link.textContent = 'Initialise units';
  initialise_units.appendChild(initialise_units_link);
  right_sidebar.appendChild(initialise_units);

  var OPF = document.createElement('p');
  OPF.style.marginLeft = '1cm';
  var OPF_link = document.createElement('a');
  OPF_link.href = 'doc_OPF.html';
  OPF_link.textContent = 'Optimal power flow';
  OPF.appendChild(OPF_link);
  right_sidebar.appendChild(OPF);

  var Dynamic = document.createElement('p');
  Dynamic.textContent = 'Dynamic';
  Dynamic.className = 'nolink_sidebar';
  right_sidebar.appendChild(Dynamic);

  var lines_and_loads = document.createElement('p');
  lines_and_loads.style.marginLeft = '1cm';
  var lines_and_loads_link = document.createElement('a');
  lines_and_loads_link.href = 'doc_lines_and_loads.html';
  lines_and_loads_link.textContent = 'Lines and loads';
  lines_and_loads.appendChild(lines_and_loads_link);
  right_sidebar.appendChild(lines_and_loads);

  var Generators = document.createElement('p');
  Generators.style.marginLeft = '1cm';
  var Generators_link = document.createElement('a');
  Generators_link.href = 'doc_Generators.html';
  Generators_link.textContent = 'Generators';
  Generators.appendChild(Generators_link);
  right_sidebar.appendChild(Generators);

  var modal_analysis = document.createElement('p');
  modal_analysis.style.marginLeft = '1cm';
  var modal_analysis_link = document.createElement('a');
  modal_analysis_link.href = 'doc_modal_analysis.html';
  modal_analysis_link.textContent = 'Modal analysis';
  modal_analysis.appendChild(modal_analysis_link);
  right_sidebar.appendChild(modal_analysis);

  // Space
  var space = document.createElement('br');
  right_sidebar.appendChild(space);

  // Add the link to VFlexP GitHub repository
  var pvflexpLink = document.createElement('p');
  var vflexpLink = document.createElement('a');
  vflexpLink.href = 'https://github.com/atomasmartin/VFlexP';
  vflexpLink.textContent = 'VFlexP GitHub repository';
  vflexpLink.style.display = 'block';
  pvflexpLink.appendChild(vflexpLink);
  right_sidebar.appendChild(pvflexpLink); 

  // Add the link to FUBM-MATPOWER GitHub repository
  var pgithubLink = document.createElement('p');
  var githubLink = document.createElement('a');
  githubLink.href = 'https://github.com/AbrahamAlvarezB/matpower-fubm';
  githubLink.textContent = 'MATPOWER-FUBM GitHub repository';
  githubLink.style.display = 'block';
  pgithubLink.appendChild(githubLink);
  right_sidebar.appendChild(pgithubLink);

  // Authors and MATLAB version
  var pAuthor1 = document.createElement('p');
  var Author1 = document.createElement('a');
  Author1.href = 'https://www.iit.comillas.edu/people/atomas';
  Author1.textContent = 'Andrés Tomás-Martín';
  pAuthor1.appendChild(Author1);
  right_sidebar.appendChild(pAuthor1);
  var pAuthor2 = document.createElement('p');
  var Author2 = document.createElement('a');
  Author2.href = 'https://www.iit.comillas.edu/people/jsuarezp';
  Author2.textContent = 'Jorge Suárez-Porras';
  pAuthor2.appendChild(Author2);
  right_sidebar.appendChild(pAuthor2);
  var pAuthor3 = document.createElement('p');
  var Author3 = document.createElement('a');
  Author3.href = 'https://www.iit.comillas.edu/people/czuluaga';
  Author3.textContent = 'Carlos-David Zuluaga-Ríos';
  pAuthor3.appendChild(Author3);
  right_sidebar.appendChild(pAuthor3);
  var pAuthor4 = document.createElement('p');
  var Author4 = document.createElement('a');
  Author4.href = 'https://www.iit.comillas.edu/people/jgaguilar';
  Author4.textContent = 'Javier García-Aguilar';
  pAuthor4.appendChild(Author4);
  right_sidebar.appendChild(pAuthor4);
  var pAuthor5 = document.createElement('p');
  var Author5 = document.createElement('a');
  Author5.href = 'https://www.iit.comillas.edu/people/aurelio';
  Author5.textContent = 'Aurelio García-Cerrada';
  pAuthor5.appendChild(Author5);
  right_sidebar.appendChild(pAuthor5);
  var pAuthor6 = document.createElement('p');
  var Author6 = document.createElement('a');
  Author6.href = 'https://www.iit.comillas.edu/people/lsigrist';
  Author6.textContent = 'Lukas Sigrist';
  pAuthor6.appendChild(Author6);
  right_sidebar.appendChild(pAuthor6);
  var pAuthor7 = document.createElement('p');
  var Author7 = document.createElement('a');
  Author7.href = 'https://www.durham.ac.uk/staff/behzad-kazemtabrizi/';
  Author7.textContent = 'Behzad Kazemtabrizi';
  pAuthor7.appendChild(Author7);
  right_sidebar.appendChild(pAuthor7);
  var pAuthor8 = document.createElement('p');
  var Author8 = document.createElement('a');
  Author8.href = 'https://www.iit.comillas.edu/people/cdepaolis';
  Author8.textContent = 'Carlo de Paolis';
  pAuthor8.appendChild(Author8);
  right_sidebar.appendChild(pAuthor8);
  var pAuthor8 = document.createElement('p');
  var Author8 = document.createElement('a');
  Author8.href = 'https://www.iit.comillas.edu/people/ravila';
  Author8.textContent = 'Régulo-Enrique Ávila-Martínez';
  pAuthor8.appendChild(Author8);
  right_sidebar.appendChild(pAuthor8);
  
  var pMATLABversion = document.createElement('p');
  var MATLABversion = document.createElement('a');
  MATLABversion.href = 'https://www.mathworks.com/products/matlab.html';
  MATLABversion.textContent = 'MATLAB R2024a';
  pMATLABversion.appendChild(MATLABversion);
  right_sidebar.appendChild(pMATLABversion);
  
  
});
