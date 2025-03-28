// Function to handle copy code button click
function copyCode(button) {
    const codeContainer = button.parentElement;
    const code = codeContainer.querySelector("pre").innerText;
    navigator.clipboard.writeText(code).then(() => {
        button.innerText = "Copied!";
        setTimeout(() => {
            button.innerText = "Copy";
        }, 2000);
    });
}

// Function to handle toggle folding
function toggleFold(button) {
    const content = button.nextElementSibling;
    if (content.classList.contains("hidden")) {
        content.classList.remove("hidden");  // Show content
        button.classList.remove("collapsed");  // Unrotate the triangle
    } else {
        content.classList.add("hidden");  // Hide content
        button.classList.add("collapsed");  // Rotate the triangle
    }
}

document.addEventListener("DOMContentLoaded", function() {
    // Add event listeners to all toggle buttons
    document.querySelectorAll('.toggle-btn').forEach(button => {
        console.log(button);  // Check if buttons are found
        button.addEventListener('click', function() {
            toggleFold(this);
        });
    });
});

// MATLAB Code Highlighter
function highlightMATLABCode() {
    const matlabKeywords = [
        'function', 'end', 'if', 'else', 'elseif', 'for', 'while', 'switch', 'case', 'otherwise', 'break',
        'continue', 'return', 'try', 'catch', 'global', 'persistent', 'set_param', 'get_param', 'disp', 'fprintf'
    ];

    const codeElements = document.querySelectorAll('.matlab-code');
    codeElements.forEach(codeElement => {
        let code = codeElement.innerText;

        // Highlight comments first
        code = code.replace(/(%[^\n]*)/g, '<span class="comment">$1</span>');

        // Highlight strings
        code = code.replace(/(['`][^'`]*['`])/g, '<span class="string">$1</span>');

        // Highlight keywords, but only outside of comments and strings
        matlabKeywords.forEach(keyword => {
            const keywordRegex = new RegExp(`\\b${keyword}\\b(?![^<]*>|[^<>]*<\\/span>)`, 'g');
            code = code.replace(keywordRegex, `<span class="keyword">${keyword}</span>`);
        });

        // Replace the content with the highlighted code
        codeElement.innerHTML = code;
    });
}

// Run the MATLAB code highlighting after page load
document.addEventListener("DOMContentLoaded", highlightMATLABCode);
