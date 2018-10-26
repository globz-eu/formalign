# from selenium import webdriver
# from formalign.settings import CHROME_DRIVER

#
# def before_scenario(context, scenario):
#     context.browser = webdriver.Chrome(CHROME_DRIVER)

def after_scenario(context, scenario):
    context.browser.quit()
